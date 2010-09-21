/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2007 H. Jasak All rights reserved
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
    Reduced Rank Extrapolation Algebraic Multigrid solver with run-time
    selection of policy and cycle

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "rreAmgSolver.H"
#include "scalarMatrices.H"
#include "DenseMatrixTools.H"
#include "FieldFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    defineTypeNameAndDebug(rreAmgSolver, 0);

    lduSolver::addsymMatrixConstructorToTable<rreAmgSolver>
        addrreAmgSolverSymMatrixConstructorToTable_;

    lduSolver::addasymMatrixConstructorToTable<rreAmgSolver>
        addrreAmgSolverAsymMatrixConstructorToTable_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from matrix and solver data stream
Foam::rreAmgSolver::rreAmgSolver
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
    amg_
    (
        matrix,
        coupleBouCoeffs,
        coupleIntCoeffs,
        interfaces,
        dict
    ),
    kDimension_(readLabel(dict.lookup("kDimension")))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::lduSolverPerformance Foam::rreAmgSolver::solve
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
    solverPerf.initialResidual() =
        gSumMag(amg_.residual(x, b, cmpt))/normFactor;

    solverPerf.finalResidual() = solverPerf.initialResidual();

    // Krylov vectors
    typedef FieldField<Field, scalar> scalarFieldField;

    scalarFieldField xSave(kDimension_ + 1);

    forAll (xSave, i)
    {
        xSave.set(i, new scalarField(x.size()));
    }

    scalarFieldField xSecondDiffs(kDimension_ - 1);

    forAll (xSecondDiffs, i)
    {
        xSecondDiffs.set(i, new scalarField(x.size()));
    }

    scalarField xRhs(x.size());

    // Matrices
    scalarSquareMatrix AN(kDimension_ - 1);
    scalarField xN(kDimension_ - 1);
    scalarField bN(kDimension_ - 1);

    // Remember last extrapolation index
    label lastExtrapolation = kDimension_;


    // Solver loop

    if (!stop(solverPerf))
    {
        scalarField psiSave(x.size());

        do
        {
            amg_.cycle(x, b, cmpt);

            label curIndex = (solverPerf.nIterations()) % (kDimension_ + 1);

            // Save x into xSave
            xSave[curIndex] = x;

            if
            (
                solverPerf.nIterations() >= lastExtrapolation
             && curIndex == kDimension_
            )
            {
                x = xSave[0];

                // Calculate first differences
                for (label i = 0; i < kDimension_; i++)
                {
                    xSave[i] = xSave[i + 1] - xSave[i];
                }

                // Calculate second differences
                for (label i = 0; i < kDimension_ - 1; i++)
                {
                    xSecondDiffs[i] = xSave[i + 1] - xSave[i];
                }

                // Perform Q-R decomposition
                DenseMatrixTools::qrDecompose
                (
                    kDimension_ - 1,
                    xSecondDiffs,
                    AN
                );

                // Solve the normal system
                for (label i = 0; i < kDimension_ - 1; i++)
                {
                    bN[i] = -gSum(xSecondDiffs[i]*xSave[0]);
                }

                DenseMatrixTools::solve(AN, xN, bN);

                // Re-asemble the solution

                for (label i = 0; i < kDimension_ - 1; i++)
                {
                    x += xN[i]*xSave[i];
                }

                // Increment last extrapolation
                lastExtrapolation = solverPerf.nIterations() + kDimension_;
            }

            // Re-calculate residual
            solverPerf.finalResidual() =
                gSumMag(amg_.residual(x, b, cmpt))/normFactor;
            solverPerf.nIterations()++;
        } while (!stop(solverPerf));
    }

    return solverPerf;
}


// ************************************************************************* //
