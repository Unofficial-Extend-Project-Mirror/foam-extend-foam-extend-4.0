/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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
    ILU solver for symmetric and asymmetric matrices.  In
    order to improve efficiency, the residual is evaluated after every
    nSweeps sweeps.

\*---------------------------------------------------------------------------*/

#include "BlockILUSolver.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from matrix and solver data stream
template<class Type>
Foam::BlockILUSolver<Type>::BlockILUSolver
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
    ilu_(matrix),
    nSweeps_(readLabel(this->dict().lookup("nSweeps")))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
typename Foam::BlockSolverPerformance<Type>
Foam::BlockILUSolver<Type>::solve
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

    Type norm = this->normFactor(x, b);

    Field<Type> residual(x.size());

    // Calculate residual
    matrix.Amul(residual, x);

    forAll (b, i)
    {
        residual[i] = b[i] - residual[i];
    }

    solverPerf.initialResidual() = cmptDivide(gSum(cmptMag(residual)),norm);
    solverPerf.finalResidual() = solverPerf.initialResidual();

    // Check convergence, solve if not converged

    if (!this->stop(solverPerf))
    {
        // Iteration loop

        Field<Type> xCorr(x.size());

        do
        {
            for (label i = 0; i < nSweeps_; i++)
            {
                ilu_.precondition(xCorr, residual);

                x += xCorr;

                solverPerf.nIterations()++;
            }

            // Re-calculate residual
            matrix.Amul(residual, x);

            forAll (b, i)
            {
                residual[i] = b[i] - residual[i];
            }

            solverPerf.finalResidual() = cmptDivide(gSum(cmptMag(residual)), norm);
        } while (!this->stop(solverPerf));
    }

    return solverPerf;
}


// ************************************************************************* //
