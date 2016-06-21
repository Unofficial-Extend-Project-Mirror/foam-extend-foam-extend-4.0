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

Class
    coarseAmgLevel

Description
    Coarse AMG level stores matrix, x and b locally

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "coarseAmgLevel.H"
#include "SubField.H"
#include "gmresSolver.H"
#include "vector2D.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::coarseAmgLevel::coarseAmgLevel
(
    autoPtr<amgMatrix> matrixPtr,
    const dictionary& dict,
    const word& policyType,
    const label groupSize,
    const label minCoarseEqns,
    const word& smootherType
)
:
    matrixPtr_(matrixPtr),
    x_(matrixPtr_->size()),
    b_(matrixPtr_->size()),
    dict_(dict),
    policyPtr_
    (
        amgPolicy::New
        (
            policyType,
            matrixPtr_->matrix(),
            groupSize,
            minCoarseEqns
        )
    ),
    smootherPtr_
    (
        lduSmoother::New
        (
            matrixPtr_->matrix(),
            matrixPtr_->coupleBouCoeffs(),
            matrixPtr_->coupleIntCoeffs(),
            matrixPtr_->interfaceFields(),
            dict
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::coarseAmgLevel::~coarseAmgLevel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalarField& Foam::coarseAmgLevel::x()
{
    return x_;
}


Foam::scalarField& Foam::coarseAmgLevel::b()
{
    return b_;
}


void Foam::coarseAmgLevel::residual
(
    const scalarField& x,
    const scalarField& b,
    const direction cmpt,
    scalarField& res
) const
{
    // Calculate residual
    matrixPtr_->matrix().Amul
    (
        res,
        x,
        matrixPtr_->coupleBouCoeffs(),
        matrixPtr_->interfaceFields(),
        cmpt
    );

    // residual = b - Ax
    forAll (b, i)
    {
        res[i] = b[i] - res[i];
    }
}


void Foam::coarseAmgLevel::restrictResidual
(
    const scalarField& x,
    const scalarField& b,
    const direction cmpt,
    scalarField& xBuffer,
    scalarField& coarseRes,
    bool preSweepsDone
) const
{
    if (preSweepsDone)
    {
        // Calculate residual
        scalarField::subField resBuf(xBuffer, x.size());

        scalarField& res = const_cast<scalarField&>
        (
            resBuf.operator const scalarField&()
        );

        residual(x, b, cmpt, res);

        policyPtr_->restrictResidual(res, coarseRes);
    }
    else
    {
        // No pre-sweeps done: x = 0 and residual = b
        policyPtr_->restrictResidual(b, coarseRes);
    }
}


void Foam::coarseAmgLevel::prolongateCorrection
(
    scalarField& x,
    const scalarField& coarseX
) const
{
    policyPtr_->prolongateCorrection(x, coarseX);
}


void Foam::coarseAmgLevel::smooth
(
    scalarField& x,
    const scalarField& b,
    const direction cmpt,
    const label nSweeps
) const
{
    smootherPtr_->smooth(x, b, cmpt, nSweeps);
}


void Foam::coarseAmgLevel::solve
(
    scalarField& x,
    const scalarField& b,
    const direction cmpt,
    const scalar tolerance,
    const scalar relTol
) const
{
    lduSolverPerformance coarseSolverPerf;

    label maxIter = Foam::min(2*policyPtr_->minCoarseEqns(), 1000);

    dictionary topLevelDict;
    topLevelDict.add("nDirections", "5");
    topLevelDict.add("minIter", 1);
    topLevelDict.add("maxIter", maxIter);
    topLevelDict.add("tolerance", tolerance);
    topLevelDict.add("relTol", relTol);

    // Avoid issues with round-off on strict tolerance setup
    // HJ, 27/Jun/2013
    x = b/matrixPtr_->matrix().diag();

    // Do not solve if the number of equations is smaller than 5
    if (policyPtr_->minCoarseEqns() < 5)
    {
        return;
    }

    // Switch off debug in top-level direct solve
    label oldDebug = lduMatrix::debug();

    if (matrixPtr_->matrix().symmetric())
    {
        topLevelDict.add("preconditioner", "Cholesky");

        coarseSolverPerf = gmresSolver
        (
            "topLevelCorr",
            matrixPtr_->matrix(),
            matrixPtr_->coupleBouCoeffs(),
            matrixPtr_->coupleIntCoeffs(),
            matrixPtr_->interfaceFields(),
            topLevelDict
        ).solve(x, b, cmpt);
    }
    else
    {
        topLevelDict.add("preconditioner", "ILU0");

        coarseSolverPerf = gmresSolver
        (
            "topLevelCorr",
            matrixPtr_->matrix(),
            matrixPtr_->coupleBouCoeffs(),
            matrixPtr_->coupleIntCoeffs(),
            matrixPtr_->interfaceFields(),
            topLevelDict
        ).solve(x, b, cmpt);
    }

    // Restore debug
    lduMatrix::debug = oldDebug;

    // Escape cases of top-level solver divergence
    if
    (
        coarseSolverPerf.nIterations() == maxIter
     && (
            coarseSolverPerf.finalResidual()
         >= coarseSolverPerf.initialResidual()
        )
    )
    {
        // Top-level solution failed.  Attempt rescue
        // HJ, 27/Jul/2013
        x = b/matrixPtr_->matrix().diag();

        // Print top level correction failure as info for user
        coarseSolverPerf.print();
    }

    if (lduMatrix::debug >= 3)
    {
        coarseSolverPerf.print();
    }
}


void Foam::coarseAmgLevel::scaleX
(
    scalarField& x,
    const scalarField& b,
    const direction cmpt,
    scalarField& xBuffer
) const
{
    // Calculate scaling
    scalarField::subField Ax(xBuffer, x.size());

    matrixPtr_->matrix().Amul
    (
        const_cast<scalarField&>(Ax.operator const scalarField&()),
        x,
        matrixPtr_->coupleBouCoeffs(),
        matrixPtr_->interfaceFields(),
        cmpt
    );

    scalar scalingFactorNum = 0.0;
    scalar scalingFactorDenom = 0.0;

    forAll(x, i)
    {
        scalingFactorNum += x[i]*b[i];
        scalingFactorDenom += x[i]*xBuffer[i]; // Note: Ax = xBuffer
    }

    vector2D scalingVector(scalingFactorNum, scalingFactorDenom);
    reduce(scalingVector, sumOp<vector2D>());

    // Scale x
    if
    (
        mag(scalingVector[0]) > GREAT
     || mag(scalingVector[1]) > GREAT
     || scalingVector[0]*scalingVector[1] <= 0
//      || mag(scalingVector[0]) < mag(scalingVector[1])
    )
    {
        // Factor = 1.0, no scaling
    }
    else
    {
        // Regular scaling with a limiter
        scalar scalingFactor =
            Foam::max
            (
                0.1,
                Foam::min
                (
                    scalingVector[0]/stabilise(scalingVector[1], SMALL),
                    10
                )
            );

        x *= scalingFactor;
    }
}


Foam::autoPtr<Foam::amgLevel> Foam::coarseAmgLevel::makeNextLevel() const
{
    if (policyPtr_->coarsen())
    {
        return autoPtr<Foam::amgLevel>
        (
            new coarseAmgLevel
            (
                policyPtr_->restrictMatrix
                (
                    matrixPtr_->coupleBouCoeffs(),
                    matrixPtr_->coupleIntCoeffs(),
                    matrixPtr_->interfaceFields()
                ),
                dict(),
                policyPtr_->type(),
                policyPtr_->groupSize(),
                policyPtr_->minCoarseEqns(),
                smootherPtr_->type()
            )
        );
    }
    else
    {
        // Final level: cannot coarsen
        return autoPtr<Foam::amgLevel>();
    }
}


// ************************************************************************* //
