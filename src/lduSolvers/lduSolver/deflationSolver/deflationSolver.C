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
    N-th order deflation solver

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "deflationSolver.H"
#include "DenseMatrixTools.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    defineTypeNameAndDebug(deflationSolver, 0);

    lduSolver::addsymMatrixConstructorToTable<deflationSolver>
        adddeflationSolverSymMatrixConstructorToTable_;

    lduSolver::addasymMatrixConstructorToTable<deflationSolver>
        adddeflationSolverAsymMatrixConstructorToTable_;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::deflationSolver::GaussSolve
(
    const label activeSize,
    const scalarSquareMatrix& A,
    scalarField& x,
    const scalarField& b
) const
{
    if (activeSize == 0)
    {
        return;
    }
    else if (activeSize == 1)
    {
        x[0] = b[0]/A[0][0];
    }
    else
    {
        // Non-trivial solution

        // Create a temporary matrix and b - they are destroyed in solver call
        scalarSquareMatrix lm(activeSize);

        for (label i = 0; i < activeSize; i++)
        {
            for (label j = 0; j < activeSize; j++)
            {
                lm[i][j] = A[i][j];
            }
        }

        scalarField lb(activeSize);

        for (label i = 0; i < activeSize; i++)
        {
            lb[i] = b[i];
        }

        DenseMatrixTools::solve(lm, x, lb);
    }
}


Foam::scalar Foam::deflationSolver::residual
(
    const scalarField& x,
    const scalarField& b,
    const direction cmpt
) const
{
    matrix_.Amul
    (
        xBuffer_,
        x,
        coupleBouCoeffs_,
        interfaces_,
        cmpt
    );

    // residual = b - Ax
    xBuffer_ -= b;

    return gSumMag(xBuffer_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from matrix and solver data stream
Foam::deflationSolver::deflationSolver
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
    rpmOrder_(readLabel(dict.lookup("rpmOrder"))),
    maxDirs_(Foam::max(readLabel(dict.lookup("maxDirections")), 2)),
    basisTol_(readScalar(dict.lookup("basisTolerance"))),
    divTol_(readScalar(dict.lookup("divergenceTolerance"))),
    nBasisSteps_(readLabel(dict.lookup("nBasisSteps"))),
    nPowerIter_(readLabel(dict.lookup("nPowerIter"))),
    xBuffer_(matrix.lduAddr().size())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::lduSolverPerformance Foam::deflationSolver::solve
(
    scalarField& x,
    const scalarField& b,
    const direction cmpt
) const
{
    // Prepare solver performance
    lduSolverPerformance solverPerf(typeName, fieldName());

    scalar normFactor = this->normFactor(x, b, cmpt);

    // Calculate initial residual
    scalar res = residual(x, b, cmpt);

    solverPerf.initialResidual() = res/normFactor;
    solverPerf.finalResidual() = solverPerf.initialResidual();

    if (!stop(solverPerf))
    {
        // Store residuals
        scalar oldRes = res;

        // Create basis vectors (Z) and deltas (Q)
        scalarFieldField Z(maxDirs_);
        scalarFieldField Q(2);
        forAll (Q, i)
        {
            Q.set(i, new scalarField(x.size()));
        }

        // Create stable (q) and unstable (p) solution
        scalarField q(x.size(), 0);
        scalarField p(x.size(), 0);

        // Field scratch space
        scalarField q0(x.size());
        scalarField t(x.size());

        // Direction size scratch space
        scalarField u(maxDirs_);
        scalarField v(maxDirs_);

        // Matrices
        scalarSquareMatrix H(maxDirs_, 0);
        scalarSquareMatrix R(maxDirs_, 0);

        // Number of directions
        label nDirs = 0;

        // Number of steps
        label nSteps = 0;

        do
        {
            nSteps++;

            // Increase the basis if diverging or converging slowly
            if
            (
                nDirs < maxDirs_
             && (
                    res > divTol_*oldRes
                 || (nSteps >= nBasisSteps_ && res > basisTol_*oldRes)
                )
            )
            {
                // Add vectors to subspace

                DenseMatrixTools::qrDecompose(2, Q, R);

                for (label qI = 0; qI < 2; qI++)
                {
                    if
                    (
                        nDirs < maxDirs_
                     && R[qI][qI] > 1e-3*Foam::max(R[0][0], R[1][1])
                    )
                    {
                        Z.set(nDirs, new scalarField(Q[qI]));
                        nDirs++;
                    }
                }

                if (lduMatrix::debug >= 2)
                {
                    Info << "Increasing the basis: " << nDirs << endl;
                }

                // Make new basis ortho-normal
                DenseMatrixTools::qrDecompose(nDirs, Z, R);

                // Do power iterations. xBuffer is used as zero rhs
                xBuffer_ = 0;

                for (label pI = 0; pI < nPowerIter_; pI++)
                {
                    for (label j = 0; j < nDirs; j++)
                    {
                        preconPtr_->precondition(Z[j], xBuffer_, cmpt);
                    }

                    // Re-orthogonalize the basis
                    DenseMatrixTools::qrDecompose(nDirs, Z, R);
                }

                // Recalculate p and q
                // u = Z^T*x
                // p = Z*u
                //   = Z*Z^T*x
                // q = x - p
                for (label j = 0; j < nDirs; j++)
                {
                    u[j] = gSumProd(Z[j], x);
                }

                p = 0;

                for (label j = 0; j < nDirs; j++)
                {
                    const scalarField& Zj = Z[j];
                    const scalar& uj = u[j];

                    for (label i = 0; i < x.size(); i++)
                    {
                        p[i] += uj*Zj[i];
                    }
                }

                // q = x - p
                subtract(q, x, p);

                // Calculate H, operating on a 2-D matrix
                // Note: xBuffer = 0 from above
                for (label j = 0; j < nDirs; j++)
                {
                    t = Z[j];

                    preconPtr_->precondition(t, xBuffer_, cmpt);

                    for (label k = 0; k < nDirs; k++)
                    {
                        if (j == k)
                        {
                            // Diagonal coefficient
                            H[j][j] = 1.0 - gSumProd(Z[k], t);
                        }
                        else
                        {
                            // Off-diagonal
                            H[k][j] = -gSumProd(Z[k], t);
                        }
                    }
                }

                // Reset steps and clear out stored differences
                forAll (Q, i)
                {
                    Q[i] = 0;
                }

                nSteps = 0;
            }

            for (label orderI = 0; orderI < rpmOrder_ + 1; orderI++)
            {
                oldRes = res;

                q0 = q;

                // Precondition the system
                t = x;

                preconPtr_->precondition(t, b, cmpt);

                // Calculate q
                // u = Z^T*t
                // q = t - Z*u
                //   = t - Z*Z^T*t
                for (label j = 0; j < nDirs; j++)
                {
                    u[j] = gSumProd(Z[j], t);
                }

                q = t;

                for (label j = 0; j < nDirs; j++)
                {
                    const scalarField& Zj = Z[j];
                    const scalar& uj = u[j];

                    for (label i = 0; i < x.size(); i++)
                    {
                        q[i] -= uj*Zj[i];
                    }
                }

                for (label i = 0; i < x.size(); i++)
                {
                    x[i] = q[i] + p[i];
                }
            }

            // Calculate p
            // u = Z^T*(t - x)
            // v = H^-1*u
            // p = p + Z*v
            //   = p + Z*H^-1*u
            //   = p + Z*H^-1*Z^T*(t - x)
            for (label j = 0; j < nDirs; j++)
            {
                // Optimisation: u[j] = sum(Z[j]*(t - x));
                scalar sum = 0;

                const scalarField& Zj = Z[j];

                for (label i = 0; i < x.size(); i++)
                {
                    sum += Zj[i]*(t[i] - x[i]);
                }

                reduce(sum, sumOp<scalar>());

                u[j] = sum;
            }

            // Solve H*v = u
            GaussSolve(nDirs, H, v, u);

            for (label j = 0; j < nDirs; j++)
            {
                const scalarField& Zj = Z[j];
                const scalar& vj = v[j];

                for (label i = 0; i < x.size(); i++)
                {
                    p[i] += vj*Zj[i];
                }
            }

            // Assemble x
            // x = p + q
            add(x, p, q);

            // Store the q-delta vector
            // Q[next] = q - q0
            subtract(Q[solverPerf.nIterations() % 2], q, q0);

            // Re-calculate residual
            res = residual(x, b, cmpt);

            solverPerf.finalResidual() = res/normFactor;
            solverPerf.nIterations()++;
        } while (!stop(solverPerf));
    }

    return solverPerf;
}


// ************************************************************************* //
