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
    Preconditioned Bi-Conjugate Gradient solver

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "coupledBicgSolver.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(coupledBicgSolver, 0);

    coupledLduSolver::addasymMatrixConstructorToTable<coupledBicgSolver>
        addBicgSolverAsymMatrixConstructorToTable_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from matrix
Foam::coupledBicgSolver::coupledBicgSolver
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

Foam::coupledSolverPerformance Foam::coupledBicgSolver::solve
(
    FieldField<Field, scalar>& x,
    const FieldField<Field, scalar>& b,
    const direction cmpt
) const
{
    // Prepare solver performance
    coupledSolverPerformance solverPerf(typeName, fieldName());

    FieldField<Field, scalar> wA(x.size());
    FieldField<Field, scalar> rA(x.size());

    forAll (x, rowI)
    {
        wA.set(rowI, new scalarField(x[rowI].size(), 0));
        rA.set(rowI, new scalarField(x[rowI].size(), 0));
    }


    // Calculate initial residual
    matrix_.Amul(wA, x, bouCoeffs_, interfaces_, cmpt);

    // Use rA as scratch space when calculating the normalisation factor
    scalar normFactor = this->normFactor(x, b, wA, rA, cmpt);

    if (coupledLduMatrix::debug >= 2)
    {
        Info<< "   Normalisation factor = " << normFactor << endl;
    }

    // Calculate residual
    // Optimised looping.  HJ, 19/Jan/2009
    forAll (rA, rowI)
    {
        const scalarField& curB = b[rowI];
        const scalarField& curWA = wA[rowI];
        scalarField& curRA = rA[rowI];

        forAll (curRA, i)
        {
            curRA[i] = curB[i] - curWA[i];
        }
    }

    solverPerf.initialResidual() = gSumMag(rA)/normFactor;
    solverPerf.finalResidual() = solverPerf.initialResidual();

    if (!solverPerf.checkConvergence(tolerance_, relTolerance_))
    {
        scalar rho = matrix_[0].great_;
        scalar rhoOld = rho;

        scalar alpha, beta, wApT;

        FieldField<Field, scalar> pA(x.size());
        FieldField<Field, scalar> pT(x.size());

        FieldField<Field, scalar> wT(x.size());
        FieldField<Field, scalar> rT = rA;

        forAll (pA, rowI)
        {
            pA.set(rowI, new scalarField(x[rowI].size(), 0));
            pT.set(rowI, new scalarField(x[rowI].size(), 0));

            wT.set(rowI, new scalarField(x[rowI].size(), 0));
        }


        do
        {
            rhoOld = rho;

            // Execute preconditioning
            preconPtr_->precondition(wA, rA, cmpt);

            // Using standard preconditioning on the transpose
            // Not sure this is correct, but the other one does not work
            // HJ, 13/Mar/2009
//             preconPtr_->preconditionT(wT, rT, cmpt);
            preconPtr_->precondition(wT, rT, cmpt);

            // Update search directions
            rho = gSumProd(wA, rT);

            beta = rho/rhoOld;

            forAll (pA, rowI)
            {
                scalarField& curPA = pA[rowI];
                const scalarField& curWA = wA[rowI];

                forAll (curPA, i)
                {
                    curPA[i] = curWA[i] + beta*curPA[i];
                }
            }

            forAll (pT, rowI)
            {
                scalarField& curPT = pT[rowI];
                const scalarField& curWT = wT[rowI];

                forAll (curPT, i)
                {
                    curPT[i] = curWT[i] + beta*curPT[i];
                }
            }

            // Update preconditioned residual
            matrix_.Amul(wA, pA, bouCoeffs_, interfaces_, cmpt);
            matrix_.Amul(wT, pT, bouCoeffs_, interfaces_, cmpt);

            wApT = gSumProd(wA, pT);


            // Check for singularity
            if (solverPerf.checkSingularity(mag(wApT)/normFactor))
            {
                break;
            }

            // Update solution and residual
            alpha = rho/wApT;

            forAll (x, rowI)
            {
                scalarField& curX = x[rowI];
                const scalarField& curPA = pA[rowI];

                forAll (curX, i)
                {
                    curX[i] += alpha*curPA[i];
                }
            }

            forAll (rA, rowI)
            {
                scalarField& curRA = rA[rowI];
                const scalarField& curWA = wA[rowI];

                forAll (curRA, i)
                {
                    curRA[i] -= alpha*curWA[i];
                }
            }

            forAll (rT, rowI)
            {
                scalarField& curRT = rT[rowI];
                const scalarField& curWT = wT[rowI];

                forAll (curRT, i)
                {
                    curRT[i] -= alpha*curWT[i];
                }
            }

            solverPerf.finalResidual() = gSumMag(rA)/normFactor;
            solverPerf.nIterations()++;
        } while (!stop(solverPerf));
    }

    return solverPerf;
}


// ************************************************************************* //
