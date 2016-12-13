/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Description
    Algebraic Multigrid solver with run-time selection of coarsening and cycle

Author
    Klas Jareteg, 2013-04-15

\*---------------------------------------------------------------------------*/

#include "BlockAMGSolver.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::BlockAMGSolver<Type>::BlockAMGSolver
(
    const word& fieldName,
    const BlockLduMatrix<Type>& matrix,
    const dictionary& dict
)
:
    BlockIterativeSolver<Type>
    (
        fieldName,
        matrix,
        dict
    ),
    amg_
    (
        matrix,
        dict
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
typename Foam::BlockSolverPerformance<Type>
Foam::BlockAMGSolver<Type>::solve
(
    Field<Type>& x,
    const Field<Type>& b
)
{
    // Prepare solver performance
    BlockSolverPerformance<Type> solverPerf
    (
        typeName,
        this->fieldName()
    );

    scalar norm = this->normFactor(x, b);

    // Calculate initial residual
    solverPerf.initialResidual() = gSum(cmptMag(amg_.residual(x, b)))/norm;
    solverPerf.finalResidual() = solverPerf.initialResidual();

    // Stop solver on divergence
    Type minResidual = solverPerf.initialResidual();
    scalar divergenceThreshold = 2;

    if (!this->stop(solverPerf))
    {
        do
        {
            amg_.cycle(x, b);

            solverPerf.finalResidual() =
                gSum(cmptMag(amg_.residual(x, b)))/norm;

            solverPerf.nIterations()++;

            // Divergence check
            if
            (
                cmptMax
                (
                    solverPerf.finalResidual()
                  - divergenceThreshold*minResidual
                ) > 0
             && solverPerf.nIterations() > 5
            )
            {
                break;
            }

            minResidual = Foam::min(minResidual, solverPerf.finalResidual());

        } while (!this->stop(solverPerf));
    }

    return solverPerf;
}


// ************************************************************************* //
