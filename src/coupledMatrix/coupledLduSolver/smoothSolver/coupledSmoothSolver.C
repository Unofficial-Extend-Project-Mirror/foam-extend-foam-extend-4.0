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
    Smoother-solver for coupled diagonal lduMatrices.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "coupledSmoothSolver.H"
#include "addToRunTimeSelectionTable.H"
#include "coupledLduSmoother.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(coupledSmoothSolver, 0);

    coupledLduSolver::addsymMatrixConstructorToTable<coupledSmoothSolver>
        addGaussSeidelSolverSymMatrixConstructorToTable_;

    coupledLduSolver::addasymMatrixConstructorToTable<coupledSmoothSolver>
        addGaussSeidelSolverAsymMatrixConstructorToTable_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from matrix
Foam::coupledSmoothSolver::coupledSmoothSolver
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
    nSweeps_(1)
{
    readControls();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::coupledSmoothSolver::readControls()
{
    coupledIterativeSolver::readControls();
    dict().readIfPresent("nSweeps", nSweeps_);
}


Foam::coupledSolverPerformance Foam::coupledSmoothSolver::solve
(
    FieldField<Field, scalar>& x,
    const FieldField<Field, scalar>& b,
    const direction cmpt
) const
{
    // Prepare solver performance
    coupledSolverPerformance solverPerf(typeName, fieldName());

    // Do a minimum number of sweeps
    // HJ, 19/Jan/2009
    if (minIter() > 0)
    {
        autoPtr<coupledLduSmoother> smootherPtr = coupledLduSmoother::New
        (
            matrix_,
            bouCoeffs_,
            intCoeffs_,
            interfaces_,
            dict()
        );

        smootherPtr->smooth
        (
            x,
            b,
            cmpt,
            minIter()
        );

        solverPerf.nIterations() += minIter();
    }

    // Now do normal sweeps.  HJ, 19/Jan/2009

    FieldField<Field, scalar> Ax(x.size());
    FieldField<Field, scalar> temp(x.size());

    forAll (x, rowI)
    {
        Ax.set(rowI, new scalarField(x[rowI].size(), 0));
        temp.set(rowI, new scalarField(x[rowI].size(), 0));
    }

    // Calculate initial residual. Note: for efficiency, swapping sign
    matrix_.Amul(Ax, x, bouCoeffs_, interfaces_, cmpt);

    scalar normFactor = this->normFactor(x, b, Ax, temp, cmpt);

    Ax -= b;

    solverPerf.initialResidual() = gSumMag(Ax)/normFactor;
    solverPerf.finalResidual() = solverPerf.initialResidual();

    if (!solverPerf.checkConvergence(tolerance_, relTolerance_))
    {
        autoPtr<coupledLduSmoother> smootherPtr =
            coupledLduSmoother::New
            (
                matrix_,
                bouCoeffs_,
                intCoeffs_,
                interfaces_,
                dict()
            );

        // Smoothing loop
        do
        {
            smootherPtr->smooth(x, b, cmpt, nSweeps_);

            // Re-calculate residual
            matrix_.Amul(Ax, x, bouCoeffs_, interfaces_, cmpt);
            Ax -= b;

            solverPerf.finalResidual() = gSumMag(Ax)/normFactor;

            solverPerf.nIterations() += nSweeps_;
        } while (!stop(solverPerf));
    }

    return solverPerf;
}


// ************************************************************************* //
