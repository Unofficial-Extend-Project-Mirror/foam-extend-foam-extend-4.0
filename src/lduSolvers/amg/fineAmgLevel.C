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
    fineAmgLevel

Description
    Finest AMG level container, using matrix and field references

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "fineAmgLevel.H"
#include "coarseAmgLevel.H"
#include "cgSolver.H"
#include "bicgStabSolver.H"

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

    dictionary fineLevelDict;
    fineLevelDict.add("minIter", 1);
    fineLevelDict.add("maxIter", 1000);
    fineLevelDict.add("tolerance", tolerance);
    fineLevelDict.add("relTol", relTol);

    // Avoid issues with round-off on strict tolerance setup
    // HJ, 27/Jun/2013
    x = b/matrix_.diag();

    if (matrix_.symmetric())
    {
        fineLevelDict.add("preconditioner", "Cholesky");

        coarseSolverPerf = cgSolver
        (
            "fineLevelCorr",
            matrix_,
            coupleBouCoeffs_,
            coupleIntCoeffs_,
            interfaceFields_,
            fineLevelDict
        ).solve(x, b, cmpt);
    }
    else
    {
        fineLevelDict.add("preconditioner", "ILU0");

        coarseSolverPerf = bicgStabSolver
        (
            "fineLevelCorr",
            matrix_,
            coupleBouCoeffs_,
            coupleIntCoeffs_,
            interfaceFields_,
            fineLevelDict
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
