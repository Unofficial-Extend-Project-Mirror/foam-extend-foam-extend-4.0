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
    Gauss-Seidel solver for symmetric and asymmetric matrices.  In
    order to improve efficiency, the residual is evaluated after every
    nSweeps sweeps.

\*---------------------------------------------------------------------------*/

#include "BlockGaussSeidelSolver.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from matrix and solver data stream
template<class Type>
Foam::BlockGaussSeidelSolver<Type>::BlockGaussSeidelSolver
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
    gs_(matrix),
    nSweeps_(readLabel(this->dict().lookup("nSweeps")))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
typename Foam::BlockSolverPerformance<Type>
Foam::BlockGaussSeidelSolver<Type>::solve
(
    Field<Type>& x,
    const Field<Type>& b
)
{
    // Create local references to avoid the spread this-> ugliness
    const BlockLduMatrix<Type>& matrix = this->matrix_;

    // Prepare solver performance
    BlockSolverPerformance<Type> solverPerf
    (
        typeName,
        this->fieldName()
    );

    scalar norm = this->normFactor(x, b);

    Field<Type> wA(x.size());

    // Calculate residual.  Note: sign of residual swapped for efficiency
    matrix.Amul(wA, x);
    wA -= b;

    solverPerf.initialResidual() = gSum(cmptMag(wA))/norm;
    solverPerf.finalResidual() = solverPerf.initialResidual();

    // Check convergence, solve if not converged

    if (!this->stop(solverPerf))
    {
        // Iteration loop

        do
        {
            for (label i = 0; i < nSweeps_; i++)
            {
                gs_.precondition(x, b);

                solverPerf.nIterations()++;
            }

            // Re-calculate residual.  Note: sign of residual swapped
            // for efficiency
            matrix.Amul(wA, x);
            wA -= b;

            solverPerf.finalResidual() = gSum(cmptMag(wA))/norm;
        } while (!this->stop(solverPerf));
    }

    return solverPerf;
}


// ************************************************************************* //
