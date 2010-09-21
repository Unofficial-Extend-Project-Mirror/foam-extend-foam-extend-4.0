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
    fineAmgLevel

Description
    Finest AMG level container, using matrix and field references

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "fineAmgLevel.H"
#include "coarseAmgLevel.H"
#include "ICCG.H"
#include "BICCG.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::fineAmgLevel::fineAmgLevel
(
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& coupleBouCoeffs,
    const FieldField<Field, scalar>& coupleIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaceFields,
    const dictionary& dict,
    const word& policyType,
    const label groupSize,
    const label minCoarseEqns,
    const word& smootherType
)
:
    matrix_(matrix),
    coupleBouCoeffs_(coupleBouCoeffs),
    coupleIntCoeffs_(coupleIntCoeffs),
    interfaceFields_(interfaceFields),
    dict_(dict),
    policyPtr_
    (
        amgPolicy::New(policyType, matrix_, groupSize, minCoarseEqns)
    ),
    smootherPtr_
    (
        lduSmoother::New
        (
            matrix,
            coupleBouCoeffs_,
            coupleIntCoeffs_,
            interfaceFields_,
            dict
        )
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalarField& Foam::fineAmgLevel::x()
{
    FatalErrorIn("scalarField& Foam::fineAmgLevel::x()")
        << "x is not available."
        << abort(FatalError);

    // Dummy return
    return const_cast<scalarField&>(scalarField::null());
}


Foam::scalarField& Foam::fineAmgLevel::b()
{
    FatalErrorIn("scalarField& Foam::fineAmgLevel::b()")
        << "b is not available."
        << abort(FatalError);

    // Dummy return
    return const_cast<scalarField&>(scalarField::null());
}


void Foam::fineAmgLevel::residual
(
    const scalarField& x,
    const scalarField& b,
    const direction cmpt,
    scalarField& res
) const
{
    matrix_.Amul
    (
        res,
        x,
        coupleBouCoeffs_,
        interfaceFields_,
        cmpt
    );

    // residual = b - Ax
    forAll (b, i)
    {
        res[i] = b[i] - res[i];
    }
}


void Foam::fineAmgLevel::restrictResidual
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
        residual(x, b, cmpt, xBuffer);
    }

    // Here x != 0.  It is assumed that the buffer will contain the residual
    // if no pre-sweeps have been done.  HJ, 4/Sep/2006
    policyPtr_->restrictResidual(xBuffer, coarseRes);
}


void Foam::fineAmgLevel::prolongateCorrection
(
    scalarField& x,
    const scalarField& coarseX
) const
{
    policyPtr_->prolongateCorrection(x, coarseX);
}


void Foam::fineAmgLevel::smooth
(
    scalarField& x,
    const scalarField& b,
    const direction cmpt,
    const label nSweeps
) const
{
    smootherPtr_->smooth(x, b, cmpt, nSweeps);
}


void Foam::fineAmgLevel::solve
(
    scalarField& x,
    const scalarField& b,
    const direction cmpt,
    const scalar tolerance,
    const scalar relTol
) const
{
    Info << "Fine level direct solution" << endl;

    lduSolverPerformance coarseSolverPerf;

    if (matrix_.symmetric())
    {
        coarseSolverPerf = ICCG
        (
            "topLevelCorr",
            matrix_,
            coupleBouCoeffs_,
            coupleIntCoeffs_,
            interfaceFields_,
            tolerance,
            relTol
        ).solve(x, b, cmpt);
    }
    else
    {
        coarseSolverPerf = BICCG
        (
            "topLevelCorr",
            matrix_,
            coupleBouCoeffs_,
            coupleIntCoeffs_,
            interfaceFields_,
            tolerance,
            relTol
        ).solve(x, b, cmpt);
    }

    if (lduMatrix::debug >= 2)
    {
        coarseSolverPerf.print();
    }
}


void Foam::fineAmgLevel::scaleX
(
    scalarField& x,
    const scalarField& b,
    const direction cmpt,
    scalarField& xBuffer
) const
{
    // Calculate scaling
    scalarField& Ax = xBuffer;

    matrix_.Amul
    (
        Ax,
        x,
        coupleBouCoeffs_,
        interfaceFields_,
        cmpt
    );

    scalar scalingFactorNum = 0.0;
    scalar scalingFactorDenom = 0.0;

    forAll(x, i)
    {
        scalingFactorNum += x[i]*b[i];
        scalingFactorDenom += x[i]*Ax[i];
    }

    vector scalingVector(scalingFactorNum, scalingFactorDenom, 0);
    reduce(scalingVector, sumOp<vector>());

    // Scale x
    if
    (
        scalingVector[0]*scalingVector[1] <= 0
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


Foam::autoPtr<Foam::amgLevel> Foam::fineAmgLevel::makeNextLevel() const
{
    if (policyPtr_->coarsen())
    {
        return autoPtr<Foam::amgLevel>
        (
            new coarseAmgLevel
            (
                policyPtr_->restrictMatrix
                (
                    coupleBouCoeffs_,
                    coupleIntCoeffs_,
                    interfaceFields_
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
