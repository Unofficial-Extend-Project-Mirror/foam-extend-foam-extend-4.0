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
    Preconditioned Conjugate Gradient solver.

\*---------------------------------------------------------------------------*/

#include "BlockCGSolver.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from matrix and solver data stream
template<class Type>
Foam::BlockCGSolver<Type>::BlockCGSolver
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
    preconPtr_
    (
        BlockLduPrecon<Type>::New
        (
            matrix,
            this->dict()
        )
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
typename Foam::BlockSolverPerformance<Type> Foam::BlockCGSolver<Type>::solve
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

    // Multiplication helper
    typename BlockCoeff<Type>::multiply mult;

    Field<Type> wA(x.size());

    // Calculate initial residual
    matrix.Amul(wA, x);
    Field<Type> rA(b - wA);

    solverPerf.initialResidual() = gSum(cmptMag(rA))/norm;
    solverPerf.finalResidual() = solverPerf.initialResidual();

    // Check convergence, solve if not converged

    if (!this->stop(solverPerf))
    {
        scalar rho = this->great_;
        scalar rhoOld = rho;

        scalar alpha, beta, wApA;

        Field<Type> pA(x.size());

        do
        {
            rhoOld = rho;

            // Execute preconditioning
            preconPtr_->precondition(wA, rA);

            // Update search directions
            rho = gSumProd(wA, rA);

            beta = rho/rhoOld;

            forAll (pA, i)
            {
                pA[i] = wA[i] + beta*pA[i];
            }

            // Update preconditioner residual
            matrix.Amul(wA, pA);

            wApA = gSumProd(wA, pA);

            // Check for singularity
            if (solverPerf.checkSingularity(mag(wApA)/norm))
            {
                break;
            }

            // Update solution and raw residual
            alpha = rho/wApA;

            forAll (x, i)
            {
                x[i] += alpha*pA[i];
            }

            forAll (rA, i)
            {
                rA[i] -= alpha*wA[i];
            }

            solverPerf.finalResidual() = gSum(cmptMag(rA))/norm;
            solverPerf.nIterations()++;
        } while (!this->stop(solverPerf));
    }

    return solverPerf;
}


// ************************************************************************* //
