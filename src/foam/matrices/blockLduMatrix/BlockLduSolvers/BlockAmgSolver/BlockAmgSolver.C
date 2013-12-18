/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
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
    Algebraic Multigrid solver with run-time selection of policy and cycle

Author
    Klas Jareteg, 2013-04-15

\*---------------------------------------------------------------------------*/

#include "BlockAmgSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from matrix and solver data stream
template<class Type>
Foam::BlockAmgSolver<Type>::BlockAmgSolver
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
Foam::BlockAmgSolver<Type>::solve
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

    // Create local references to avoid the spread this-> ugliness
    const BlockLduMatrix<Type>& matrix = this->matrix_;

    scalar norm = this->normFactor(x, b);

    Field<Type> wA(x.size());

    // Calculate residual.  Note: sign of residual swapped for efficiency
    matrix.Amul(wA, x);
    wA -= b;

    solverPerf.initialResidual() = gSum(cmptMag(wA))/norm;
    solverPerf.finalResidual() = solverPerf.initialResidual();

    if (!solverPerf.checkConvergence(this->tolerance(), this->relTolerance()))
    {
        do
        {
            amg_.cycle(x, b);

            // Re-calculate residual.  Note: sign of residual swapped
            // for efficiency
            matrix.Amul(wA, x);
            wA -= b;

            solverPerf.finalResidual() = gSum(cmptMag(wA))/norm;
            solverPerf.nIterations()++;

        } while (!this->stop(solverPerf));
    }

    return solverPerf;
}


// ************************************************************************* //
