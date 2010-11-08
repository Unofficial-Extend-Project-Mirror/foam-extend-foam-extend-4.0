/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-6 H. Jasak All rights reserved
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Description
    Preconditioned Generalised Minimal Residual solver with
    run-time selectable preconditioning

\*---------------------------------------------------------------------------*/

#include "BlockGMRESSolver.H"
#include "FieldField.H"
#include "scalarMatrices.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::BlockGMRESSolver<Type>::givensRotation
(
    const scalar& h,
    const scalar& beta,
    scalar& c,
    scalar& s
) const
{
    if (beta == 0)
    {
        c = 1;
        s = 0;
    }
    else if (mag(beta) > mag(h))
    {
        scalar tau = -h/beta;
        s = 1.0/Foam::sqrt(1.0 + sqr(tau));
        c = s*tau;
    }
    else
    {
        scalar tau = -beta/h;
        c = 1.0/Foam::sqrt(1.0 + sqr(tau));
        s = c*tau;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from matrix and solver data stream
template<class Type>
Foam::BlockGMRESSolver<Type>::BlockGMRESSolver
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
    ),
    nDirs_(readLabel(this->dict().lookup("nDirections")))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
typename Foam::BlockSolverPerformance<Type>
Foam::BlockGMRESSolver<Type>::solve
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

    if (!solverPerf.checkConvergence(this->tolerance(), this->relTolerance()))
    {
        // Create the Hesenberg matrix
        scalarSquareMatrix H(nDirs_, 0);

        // Create y and b for Hessenberg matrix
        scalarField yh(nDirs_, 0);
        scalarField bh(nDirs_ + 1, 0);

        // Givens rotation vectors
        scalarField c(nDirs_, 0);
        scalarField s(nDirs_, 0);

        // Allocate Krylov space vectors
        FieldField<Field, Type> V(nDirs_ + 1);

        forAll (V, i)
        {
            V.set(i, new Field<Type>(x.size(), pTraits<Type>::zero));
        }

        do
        {
            // Execute preconditioning
            preconPtr_->precondition(wA, rA);

            // Calculate beta and scale first vector
            scalar beta = Foam::sqrt(gSumProd(wA, wA));

            // Set initial rhs and bh[0] = beta
            bh = 0;
            bh[0] = beta;

            for (label i = 0; i < nDirs_; i++)
            {
                // Set search direction
                V[i] = wA;
                V[i] /= beta;

                // Arnoldi's method
                matrix.Amul(rA, V[i]);

                // Execute preconditioning
                preconPtr_->precondition(wA, rA);

                for (label j = 0; j <= i; j++)
                {
                    beta = gSumProd(wA, V[j]);

                    H[j][i] = beta;

                    forAll (wA, wI)
                    {
                        wA[wI] -= beta*V[j][wI];
                    }
                }

                beta = Foam::sqrt(gSumProd(wA, wA));

                // Apply previous Givens rotations to new column of H.
                for (label j = 0; j < i; j++)
                {
                    const scalar Hji = H[j][i];
                    H[j][i] = c[j]*Hji - s[j]*H[j + 1][i];
                    H[j + 1][i] = s[j]*Hji + c[j]*H[j + 1][i];
                }

                // Apply Givens rotation to current row.
                givensRotation(H[i][i], beta, c[i], s[i]);

                const scalar bhi = bh[i];
                bh[i] = c[i]*bhi - s[i]*bh[i + 1];
                bh[i + 1] = s[i]*bhi + c[i]*bh[i + 1];
                H[i][i] = c[i]*H[i][i] - s[i]*beta;
            }

            // Back substitute to solve Hy = b
            for (label i = nDirs_ - 1; i >= 0; i--)
            {
                scalar sum = bh[i];

                for (label j = i + 1; j < nDirs_; j++)
                {
                    sum -= H[i][j]*yh[j];
                }

                yh[i] = sum/H[i][i];
            }

            // Update solution

            for (label i = 0; i < nDirs_; i++)
            {
                const Field<Type>& Vi = V[i];
                const scalar& yi = yh[i];

                forAll (x, xI)
                {
                    x[xI] += yi*Vi[xI];
                }
            }

            // Re-calculate the residual
            matrix.Amul(wA, x);

            forAll (rA, raI)
            {
                rA[raI] = b[raI] - wA[raI];
            }

            solverPerf.finalResidual() = gSum(cmptMag(rA))/norm;
            solverPerf.nIterations()++;
        } while (!stop(solverPerf));
    }

    return solverPerf;
}


// ************************************************************************* //
