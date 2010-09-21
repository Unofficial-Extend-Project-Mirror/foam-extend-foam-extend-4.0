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
    symGaussSeidelPrecon

Description
    Symmetric Gauss-Seidel preconditioning

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "symGaussSeidelPrecon.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(symGaussSeidelPrecon, 0);

    lduPreconditioner::
        addsymMatrixConstructorToTable<symGaussSeidelPrecon>
        addsymGaussSeidelPreconditionerSymMatrixConstructorToTable_;

    lduPreconditioner::
        addasymMatrixConstructorToTable<symGaussSeidelPrecon>
        addsymGaussSeidelPreconditionerAsymMatrixConstructorToTable_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::symGaussSeidelPrecon::symGaussSeidelPrecon
(
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& coupleBouCoeffs,
    const FieldField<Field, scalar>& coupleIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces
)
:
    lduPreconditioner
    (
        matrix,
        coupleBouCoeffs,
        coupleIntCoeffs,
        interfaces
    ),
    mBouCoeffs_(coupleBouCoeffs.size()),
    bPrime_(matrix.lduAddr().size())
{
    forAll(mBouCoeffs_, i)
    {
        if (interfaces_.set(i))
        {
            mBouCoeffs_.set(i, -coupleBouCoeffs_[i]);
        }
    }
}


Foam::symGaussSeidelPrecon::symGaussSeidelPrecon
(
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& coupleBouCoeffs,
    const FieldField<Field, scalar>& coupleIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const dictionary& dict
)
:
    lduPreconditioner
    (
        matrix,
        coupleBouCoeffs,
        coupleIntCoeffs,
        interfaces
    ),
    mBouCoeffs_(coupleBouCoeffs.size()),
    bPrime_(matrix.lduAddr().size())
{
    forAll(mBouCoeffs_, i)
    {
        if (interfaces_.set(i))
        {
            mBouCoeffs_.set(i, -coupleBouCoeffs_[i]);
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::symGaussSeidelPrecon::precondition
(
    scalarField& x,
    const scalarField& b,
    const direction cmpt
) const
{
    // Execute preconditioning
    if (matrix_.diagonal())
    {
        x = b/matrix_.diag();
    }
    else if (matrix_.symmetric() || matrix_.asymmetric())
    {
        scalar* __restrict__ xPtr = x.begin();

        const scalar* __restrict__ diagPtr = matrix_.diag().begin();

        scalar* __restrict__ bPrimePtr = bPrime_.begin();

        const label* const __restrict__ uPtr =
            matrix_.lduAddr().upperAddr().begin();

        const label* const __restrict__ ownStartPtr =
            matrix_.lduAddr().ownerStartAddr().begin();

        const scalar* const __restrict__ lowerPtr = matrix_.lower().begin();
        const scalar* const __restrict__ upperPtr = matrix_.upper().begin();

        const label nRows = x.size();

        label fStart, fEnd;

        bPrime_ = b;

        // Coupled boundary update
        {
            matrix_.initMatrixInterfaces
            (
                mBouCoeffs_,
                interfaces_,
                x,
                bPrime_,
                cmpt
            );

            matrix_.updateMatrixInterfaces
            (
                mBouCoeffs_,
                interfaces_,
                x,
                bPrime_,
                cmpt
            );
        }

        // Forward sweep
        for (register label rowI = 0; rowI < nRows; rowI++)
        {
            // lRow is equal to rowI
            scalar& curX = xPtr[rowI];

            // Grab the accumulated neighbour side
            curX = bPrimePtr[rowI];

            // Start and end of this row
            fStart = ownStartPtr[rowI];
            fEnd = ownStartPtr[rowI + 1];

            // Accumulate the owner product side
            for (register label curCoeff = fStart; curCoeff < fEnd; curCoeff++)
            {
                curX -= upperPtr[curCoeff]*xPtr[uPtr[curCoeff]];
            }

            // Finish current x
            curX /= diagPtr[rowI];

            // Distribute the neighbour side using current x
            for (register label curCoeff = fStart; curCoeff < fEnd; curCoeff++)
            {
                bPrimePtr[uPtr[curCoeff]] -= lowerPtr[curCoeff]*curX;
            }
        }

        // Reverse sweep
        for (register label rowI = nRows - 1; rowI >= 0; rowI--)
        {
            // lRow is equal to rowI
            scalar& curX = xPtr[rowI];

            // Grab the accumulated neighbour side
            curX = bPrimePtr[rowI];

            // Start and end of this row
            fStart = ownStartPtr[rowI];
            fEnd = ownStartPtr[rowI + 1];

            // Accumulate the owner product side
            for (register label curCoeff = fStart; curCoeff < fEnd; curCoeff++)
            {
                curX -= upperPtr[curCoeff]*xPtr[uPtr[curCoeff]];
            }

            // Finish current x
            curX /= diagPtr[rowI];

            // Distribute the neighbour side using current x
            for (register label curCoeff = fStart; curCoeff < fEnd; curCoeff++)
            {
                bPrimePtr[uPtr[curCoeff]] -= lowerPtr[curCoeff]*curX;
            }
        }
    }
}


void Foam::symGaussSeidelPrecon::preconditionT
(
    scalarField& x,
    const scalarField& b,
    const direction cmpt
) const
{
    // Execute preconditioning
    if (matrix_.diagonal())
    {
        x = b/matrix_.diag();
    }
    else if (matrix_.symmetric() || matrix_.asymmetric())
    {
        scalar* __restrict__ xPtr = x.begin();

        const scalar* __restrict__ diagPtr = matrix_.diag().begin();

        scalar* __restrict__ bPrimePtr = bPrime_.begin();

        const label* const __restrict__ uPtr =
            matrix_.lduAddr().upperAddr().begin();

        const label* const __restrict__ ownStartPtr =
            matrix_.lduAddr().ownerStartAddr().begin();

        const scalar* const __restrict__ lowerPtr = matrix_.lower().begin();
        const scalar* const __restrict__ upperPtr = matrix_.upper().begin();

        const label nRows = x.size();

        label fStart, fEnd;

        bPrime_ = b;

        // Coupled boundary update
        {
            matrix_.initMatrixInterfaces
            (
                mBouCoeffs_,
                interfaces_,
                x,
                bPrime_,
                cmpt
            );

            matrix_.updateMatrixInterfaces
            (
                mBouCoeffs_,
                interfaces_,
                x,
                bPrime_,
                cmpt
            );
        }

        // Forward sweep
        for (register label rowI = 0; rowI < nRows; rowI++)
        {
            // lRow is equal to rowI
            scalar& curX = xPtr[rowI];

            // Grab the accumulated neighbour side
            curX = bPrimePtr[rowI];

            // Start and end of this row
            fStart = ownStartPtr[rowI];
            fEnd = ownStartPtr[rowI + 1];

            // Accumulate the owner product side
            for (register label curCoeff = fStart; curCoeff < fEnd; curCoeff++)
            {
                // Transpose multiplication.  HJ, 10/Jul/2007
                curX -= lowerPtr[curCoeff]*xPtr[uPtr[curCoeff]];
            }

            // Finish current x
            curX /= diagPtr[rowI];

            // Distribute the neighbour side using current x
            for (register label curCoeff = fStart; curCoeff < fEnd; curCoeff++)
            {
                // Transpose multiplication.  HJ, 10/Jul/2007
                bPrimePtr[uPtr[curCoeff]] -= upperPtr[curCoeff]*curX;
            }
        }

        // Reverse sweep
        for (register label rowI = nRows - 1; rowI >= 0; rowI--)
        {
            // lRow is equal to rowI
            scalar& curX = xPtr[rowI];

            // Grab the accumulated neighbour side
            curX = bPrimePtr[rowI];

            // Start and end of this row
            fStart = ownStartPtr[rowI];
            fEnd = ownStartPtr[rowI + 1];

            // Accumulate the owner product side
            for (register label curCoeff = fStart; curCoeff < fEnd; curCoeff++)
            {
                // Transpose multiplication.  HJ, 10/Jul/2007
                curX -= lowerPtr[curCoeff]*xPtr[uPtr[curCoeff]];
            }

            // Finish current x
            curX /= diagPtr[rowI];

            // Distribute the neighbour side using current x
            for (register label curCoeff = fStart; curCoeff < fEnd; curCoeff++)
            {
                // Transpose multiplication.  HJ, 10/Jul/2007
                bPrimePtr[uPtr[curCoeff]] -= upperPtr[curCoeff]*curX;
            }
        }
    }
}


// ************************************************************************* //
