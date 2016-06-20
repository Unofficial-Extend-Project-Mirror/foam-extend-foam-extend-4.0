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
    Minimal Polynomial Extrapolation Algebraic Multigrid solver with run-time
    selection of policy and cycle

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "mpeAmgSolver.H"
#include "scalarMatrices.H"
#include "DenseMatrixTools.H"
#include "FieldFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    defineTypeNameAndDebug(mpeAmgSolver, 0);

    lduSolver::addsymMatrixConstructorToTable<mpeAmgSolver>
        addmpeAmgSolverSymMatrixConstructorToTable_;

    lduSolver::addasymMatrixConstructorToTable<mpeAmgSolver>
        addmpeAmgSolverAsymMatrixConstructorToTable_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from matrix and solver data stream
Foam::mpeAmgSolver::mpeAmgSolver
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

Foam::lduSolverPerformance Foam::mpeAmgSolver::solve
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


    // Solver loop

    if (!stop(solverPerf))
    {
        // Krylov vectors
        typedef FieldField<Field, scalar> scalarFieldField;

        scalarFieldField xSave(kDimension_ + 1);

        forAll (xSave, i)
        {
            xSave.set(i, new scalarField(x.size()));
        }

        scalarFieldField xFirstDiffs(kDimension_ - 1);

        forAll (xFirstDiffs, i)
        {
            xFirstDiffs.set(i, new scalarField(x.size()));
        }

        scalarField xRhs(x.size());

        // Matrices
        scalarSquareMatrix AN(kDimension_ - 1);
        scalarField xN(kDimension_ - 1);
        scalarField bN(kDimension_ - 1);

        // Remember last extrapolation index
        label lastExtrapolation = kDimension_;

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
                // Calculate differences
                for (label i = 0; i < kDimension_ - 1; i++)
                {
                    xFirstDiffs[i] = xSave[i] - xSave[kDimension_];
                }

                // Grab right-hand side
                xRhs = xSave[kDimension_] - xSave[kDimension_ - 1];

                // Perform Q-R decomposition
                DenseMatrixTools::qrDecompose
                (
                    kDimension_ - 1,
                    xFirstDiffs,
                    AN
                );

                // Solve the normal system
                for (label i = 0; i < kDimension_ - 1; i++)
                {
                    bN[i] = -gSum(xFirstDiffs[i]*xRhs);
                }

                DenseMatrixTools::solve(AN, xN, bN);

                // Reset last interpolation factor
                xN[kDimension_ - 2] = 1.0;

                // Re-normalize
                xN /= sum(xN);


                // Re-asemble the solution.  HJ, variants?
                x = xSave[kDimension_];

                for (label i = 0; i < kDimension_ - 1; i++)
                {
                    xFirstDiffs[0] = xSave[i + 1] - xSave[i];
                    x += xN[i]*xFirstDiffs[0];
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
