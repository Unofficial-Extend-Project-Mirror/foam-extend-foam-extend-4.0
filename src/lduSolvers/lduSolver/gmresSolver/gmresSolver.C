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
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    Preconditioned Generalised Minimal Residual solver with
    run-time selectable preconditioning

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "gmresSolver.H"
#include "scalarMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(gmresSolver, 0);

    lduSolver::addsymMatrixConstructorToTable<gmresSolver>
        addgmresSolverSymMatrixConstructorToTable_;

    lduSolver::addasymMatrixConstructorToTable<gmresSolver>
        addgmresSolverAsymMatrixConstructorToTable_;

}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::gmresSolver::givensRotation
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
Foam::gmresSolver::gmresSolver
(
    const word& fieldName,
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& coupleBouCoeffs,
    const FieldField<Field, scalar>& coupleIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const dictionary& dict
)
:
    lduSolver
    (
        fieldName,
        matrix,
        coupleBouCoeffs,
        coupleIntCoeffs,
        interfaces,
        dict
    ),
    preconPtr_
    (
        lduPreconditioner::New
        (
            matrix,
            coupleBouCoeffs,
            coupleIntCoeffs,
            interfaces,
            dict
        )
    ),
    nDirs_(readLabel(dict.lookup("nDirections")))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::lduSolverPerformance Foam::gmresSolver::solve
(
    scalarField& x,
    const scalarField& b,
    const direction cmpt
) const
{
    // Prepare solver performance
    lduSolverPerformance solverPerf(typeName, fieldName());

    scalarField wA(x.size());
    scalarField rA(x.size());

    // Calculate initial residual
    matrix_.Amul(wA, x, coupleBouCoeffs_, interfaces_, cmpt);

    // Use rA as scratch space when calculating the normalisation factor
    scalar normFactor = this->normFactor(x, b, wA, rA, cmpt);

    if (lduMatrix::debug >= 2)
    {
        Info<< "   Normalisation factor = " << normFactor << endl;
    }

    // Calculate residual
    forAll (rA, i)
    {
        rA[i] = b[i] - wA[i];
    }

    solverPerf.initialResidual() = gSumMag(rA)/normFactor;
    solverPerf.finalResidual() = solverPerf.initialResidual();

    if (!stop(solverPerf))
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
        FieldField<Field, scalar> V(nDirs_ + 1);

        forAll (V, i)
        {
            V.set(i, new scalarField(x.size(), 0));
        }

        do
        {
            // Execute preconditioning
            preconPtr_->precondition(wA, rA, cmpt);

            // Calculate beta and scale first vector
            scalar beta = Foam::sqrt(gSumSqr(wA));

            // Set initial rhs and bh[0] = beta
            bh = 0;
            bh[0] = beta;

            for (label i = 0; i < nDirs_; i++)
            {
                // Set search direction
                V[i] = wA;
                V[i] /= beta;

                // Arnoldi's method
                matrix_.Amul(rA, V[i], coupleBouCoeffs_, interfaces_, cmpt);

                // Execute preconditioning
                preconPtr_->precondition(wA, rA, cmpt);

                for (label j = 0; j <= i; j++)
                {
                    beta = gSumProd(wA, V[j]);

                    H[j][i] = beta;

                    forAll (wA, wI)
                    {
                        wA[wI] -= beta*V[j][wI];
                    }
                }

                beta = Foam::sqrt(gSumSqr(wA));

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
                const scalarField& Vi = V[i];
                const scalar& yi = yh[i];

                forAll (x, psiI)
                {
                    x[psiI] += yi*Vi[psiI];
                }
            }

            // Re-calculate the residual
            matrix_.Amul(wA, x, coupleBouCoeffs_, interfaces_, cmpt);

            forAll (rA, raI)
            {
                rA[raI] = b[raI] - wA[raI];
            }

            solverPerf.finalResidual() = gSumMag(rA)/normFactor;
            solverPerf.nIterations()++;
        } while (!stop(solverPerf));
    }

    return solverPerf;
}


// ************************************************************************* //
