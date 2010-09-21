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

Class
    coarseAmgLevel

Description
    Coarse AMG level stores matrix, x and b locally

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "coarseAmgLevel.H"
#include "SubField.H"
#include "ICCG.H"
#include "BICCG.H"
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

        scalarField& res = reinterpret_cast<scalarField&>(resBuf);

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

    if (matrixPtr_->matrix().symmetric())
    {
        coarseSolverPerf = ICCG
        (
            "topLevelCorr",
            matrixPtr_->matrix(),
            matrixPtr_->coupleBouCoeffs(),
            matrixPtr_->coupleIntCoeffs(),
            matrixPtr_->interfaceFields(),
            tolerance,
            relTol
        ).solve(x, b, cmpt);
    }
    else
    {
        coarseSolverPerf = BICCG
        (
            "topLevelCorr",
            matrixPtr_->matrix(),
            matrixPtr_->coupleBouCoeffs(),
            matrixPtr_->coupleIntCoeffs(),
            matrixPtr_->interfaceFields(),
            tolerance,
            relTol
        ).solve(x, b, cmpt);
    }

    if (lduMatrix::debug >= 2)
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
        reinterpret_cast<scalarField&>(Ax),
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
     || mag(scalingVector[0]) < mag(scalingVector[1])
    )
    {
        // Factor = 1.0, no scaling
    }
    else if (mag(scalingVector[0]) > 2*mag(scalingVector[1]))
    {
        // Max factor = 2
        x *= 2.0;
    }
    else
    {
        // Regular scaling
        x *= scalingVector[0]/stabilise(scalingVector[1], SMALL);
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
