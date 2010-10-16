/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
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
    Preconditioned Stabilised Bi-Conjugate Gradient solver

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "coupledBicgStabSolver.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(coupledBicgStabSolver, 0);

    coupledLduSolver::addsymMatrixConstructorToTable<coupledBicgStabSolver>
        addBicgStabSolverSymMatrixConstructorToTable_;
    coupledLduSolver::addasymMatrixConstructorToTable<coupledBicgStabSolver>
        addBicgStabSolverAsymMatrixConstructorToTable_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from matrix
Foam::coupledBicgStabSolver::coupledBicgStabSolver
(
    const word& fieldName,
    const coupledLduMatrix& matrix,
    const PtrList<FieldField<Field, scalar> >& bouCoeffs,
    const PtrList<FieldField<Field, scalar> >& intCoeffs,
    const lduInterfaceFieldPtrsListList& interfaces,
    const dictionary& solverData
)
:
    coupledIterativeSolver
    (
        fieldName,
        matrix,
        bouCoeffs,
        intCoeffs,
        interfaces,
        solverData
    ),
    preconPtr_
    (
        coupledLduPrecon::New
        (
            matrix,
            bouCoeffs,
            intCoeffs,
            interfaces,
            dict()
        )
    )
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::coupledSolverPerformance Foam::coupledBicgStabSolver::solve
(
    FieldField<Field, scalar>& x,
    const FieldField<Field, scalar>& b,
    const direction cmpt
) const
{
    // Prepare solver performance
    coupledSolverPerformance solverPerf(typeName, fieldName());

    FieldField<Field, scalar> p(x.size());
    FieldField<Field, scalar> r(x.size());

    forAll (x, rowI)
    {
        p.set(rowI, new scalarField(x[rowI].size(), 0));
        r.set(rowI, new scalarField(x[rowI].size(), 0));
    }


    // Calculate initial residual
    matrix_.Amul(p, x, bouCoeffs_, interfaces_, cmpt);

    // Use r as scratch space when calculating the normalisation factor
    scalar normFactor = this->normFactor(x, b, p, r, cmpt);

    if (coupledLduMatrix::debug >= 2)
    {
        Info<< "   Normalisation factor = " << normFactor << endl;
    }

    // Calculate residual
    forAll (r, rowI)
    {
        const scalarField& curB = b[rowI];
        const scalarField& curP = p[rowI];
        scalarField& curR = r[rowI];

        forAll (curR, i)
        {
            curR[i] = curB[i] - curP[i];
        }
    }

    solverPerf.initialResidual() = gSumMag(r)/normFactor;
    solverPerf.finalResidual() = solverPerf.initialResidual();

    if (!solverPerf.checkConvergence(tolerance_, relTolerance_))
    {
        scalar rho = matrix_[0].great_;
        scalar rhoOld = rho;

        scalar alpha = 0;
        scalar omega = matrix_[0].great_;
        scalar beta;

        p = 0;
        FieldField<Field, scalar> ph(x.size());
        FieldField<Field, scalar> v(x.size());
        FieldField<Field, scalar> s(x.size());
        FieldField<Field, scalar> sh(x.size());
        FieldField<Field, scalar> t(x.size());

        forAll (ph, rowI)
        {
            ph.set(rowI, new scalarField(x[rowI].size(), 0));
            v.set(rowI, new scalarField(x[rowI].size(), 0));
            s.set(rowI, new scalarField(x[rowI].size(), 0));
            sh.set(rowI, new scalarField(x[rowI].size(), 0));
            t.set(rowI, new scalarField(x[rowI].size(), 0));
        }

        // Calculate transpose residual
        FieldField<Field, scalar> rw(r);

        do
        {
            rhoOld = rho;

            // Update search directions
            rho = gSumProd(rw, r);

            beta = rho/rhoOld*(alpha/omega);

            // Restart if breakdown occurs
            if (rho == 0)
            {
                rw = r;
                rho = gSumProd(rw, r);

                alpha = 0;
                omega = 0;
                beta = 0;
            }

            forAll (p, rowI)
            {
                scalarField& curP = p[rowI];
                const scalarField& curR = r[rowI];
                const scalarField& curV = v[rowI];

                forAll (curP, i)
                {
                    curP[i] = curR[i] + beta*curP[i] - beta*omega*curV[i];
                }
            }

            // Execute preconditioning
            preconPtr_->precondition(ph, p, cmpt);
            matrix_.Amul(v, ph, bouCoeffs_, interfaces_, cmpt);
            alpha = rho/gSumProd(rw, v);

            forAll (s, rowI)
            {
                scalarField& curS = s[rowI];
                const scalarField& curR = r[rowI];
                const scalarField& curV = v[rowI];

                forAll (curS, i)
                {
                    curS[i] = curR[i] - alpha*curV[i];
                }
            }

            // Execute preconditioning transpose
            preconPtr_->preconditionT(sh, s, cmpt);
            matrix_.Amul(t, sh, bouCoeffs_, interfaces_, cmpt);
            omega = gSumProd(t, s)/gSumProd(t, t);

            // Update solution and residual
            forAll (x, rowI)
            {
                scalarField& curX = x[rowI];
                const scalarField& curPh = ph[rowI];
                const scalarField& curSh = sh[rowI];

                forAll (curX, i)
                {
                    curX[i] = curX[i] + alpha*curPh[i] + omega*curSh[i];
                }
            }

            forAll (r, rowI)
            {
                scalarField& curR = r[rowI];
                const scalarField& curS = s[rowI];
                const scalarField& curT = t[rowI];

                forAll (curR, i)
                {
                    curR[i] = curS[i] - omega*curT[i];
                }
            }

            solverPerf.finalResidual() = gSumMag(r)/normFactor;
            solverPerf.nIterations()++;
        } while (!stop(solverPerf));
    }

    return solverPerf;
}


// ************************************************************************* //
