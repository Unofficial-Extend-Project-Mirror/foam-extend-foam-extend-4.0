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
    Preconditioned Conjugate Gradient solver

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "coupledCgSolver.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(coupledCgSolver, 0);

    coupledLduSolver::addsymMatrixConstructorToTable<coupledCgSolver>
        addCgSolverSymMatrixConstructorToTable_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from matrix
Foam::coupledCgSolver::coupledCgSolver
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

Foam::coupledSolverPerformance Foam::coupledCgSolver::solve
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

    // Optimised looping.  HJ, 19/Jan/2009
    forAll (rA, i)
    {
        const scalarField& bi = b[i];
        const scalarField& wAi = wA[i];
        scalarField& rAi = rA[i];

        forAll (rAi, j)
        {
            rAi[j] = bi[j] - wAi[j];
        }
    }

    solverPerf.initialResidual() = gSumMag(rA)/normFactor;
    solverPerf.finalResidual() = solverPerf.initialResidual();

    if (!solverPerf.checkConvergence(tolerance_, relTolerance_))
    {
        scalar rho = matrix_[0].great_;
        scalar rhoOld = rho;

        scalar alpha, beta, wApA;

        FieldField<Field, scalar> pA(x.size());

        forAll (pA, rowI)
        {
            pA.set(rowI, new scalarField(x[rowI].size(), 0));
        }

        do
        {
            rhoOld = rho;

            // Execute preconditioning
            preconPtr_->precondition(wA, rA, cmpt);

            // Update search directions
            rho = gSumProd(wA, rA);

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

            // Update preconditioned residual
            matrix_.Amul(wA, pA, bouCoeffs_, interfaces_, cmpt);

            wApA = gSumProd(wA, pA);


            // Check for singularity
            if (solverPerf.checkSingularity(mag(wApA)/normFactor))
            {
                break;
            }

            // Update solution and residual
            alpha = rho/wApA;

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

            solverPerf.finalResidual() = gSumMag(rA)/normFactor;
            solverPerf.nIterations()++;
        } while (!stop(solverPerf));
    }

    return solverPerf;
}


// ************************************************************************* //
