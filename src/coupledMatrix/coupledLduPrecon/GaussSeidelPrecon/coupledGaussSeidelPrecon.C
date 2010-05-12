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

Class
    coupledGaussSeidelPrecon

Description
    GaussSeidel preconditioning

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

\*----------------------------------------------------------------------------*/

#include "coupledGaussSeidelPrecon.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(coupledGaussSeidelPrecon, 0);

    addToRunTimeSelectionTable
    (
        coupledLduPrecon,
        coupledGaussSeidelPrecon,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::coupledGaussSeidelPrecon::forwardSweep
(
    const lduMatrix& matrix,
    scalarField& x,
    scalarField& bPrime
) const
{
    const scalarField& diag = matrix.diag();
    const scalarField& lower = matrix.lower();
    const scalarField& upper = matrix.upper();

    const labelList& upperAddr = matrix.lduAddr().upperAddr();
    const labelList& ownStartAddr = matrix.lduAddr().ownerStartAddr();

    const label nRows = x.size();
    label fStart, fEnd;

    for (register label rowI = 0; rowI < nRows; rowI++)
    {
        // lRow is equal to rowI
        scalar& curX = x[rowI];

        // Grab the accumulated neighbour side
        curX = bPrime[rowI];

        // Start and end of this row
        fStart = ownStartAddr[rowI];
        fEnd = ownStartAddr[rowI + 1];

        // Accumulate the owner product side
        for (register label curCoeff = fStart; curCoeff < fEnd; curCoeff++)
        {
            curX -= upper[curCoeff]*x[upperAddr[curCoeff]];
        }

        // Finish current x
        curX /= diag[rowI];

        // Distribute the neighbour side using current x
        for (register label curCoeff = fStart; curCoeff < fEnd; curCoeff++)
        {
            bPrime[upperAddr[curCoeff]] -= lower[curCoeff]*curX;
        }
    }
}


void Foam::coupledGaussSeidelPrecon::reverseSweep
(
    const lduMatrix& matrix,
    scalarField& x,
    scalarField& bPrime
) const
{
    const scalarField& diag = matrix.diag();
    const scalarField& lower = matrix.lower();
    const scalarField& upper = matrix.upper();

    const labelList& upperAddr = matrix.lduAddr().upperAddr();
    const labelList& ownStartAddr = matrix.lduAddr().ownerStartAddr();

    const label nRows = x.size();
    label fStart, fEnd;

    for (register label rowI = nRows - 1; rowI >= 0; rowI--)
    {
        // lRow is equal to rowI
        scalar& curX = x[rowI];

        // Grab the accumulated neighbour side
        curX = bPrime[rowI];

        // Start and end of this row
        fStart = ownStartAddr[rowI];
        fEnd = ownStartAddr[rowI + 1];

        // Accumulate the owner product side
        for (register label curCoeff = fStart; curCoeff < fEnd; curCoeff++)
        {
            curX -= upper[curCoeff]*x[upperAddr[curCoeff]];
        }

        // Finish current x
        curX /= diag[rowI];

        // Distribute the neighbour side using current x
        for (register label curCoeff = fStart; curCoeff < fEnd; curCoeff++)
        {
            bPrime[upperAddr[curCoeff]] -= lower[curCoeff]*curX;
        }
    }
}


void Foam::coupledGaussSeidelPrecon::forwardSweepTranspose
(
    const lduMatrix& matrix,
    scalarField& x,
    scalarField& bPrime
) const
{
    const scalarField& diag = matrix.diag();
    const scalarField& lower = matrix.lower();
    const scalarField& upper = matrix.upper();

    const labelList& upperAddr = matrix.lduAddr().upperAddr();
    const labelList& ownStartAddr = matrix.lduAddr().ownerStartAddr();

    const label nRows = x.size();
    label fStart, fEnd;

    for (register label rowI = 0; rowI < nRows; rowI++)
    {
        // lRow is equal to rowI
        scalar& curX = x[rowI];

        // Grab the accumulated neighbour side
        curX = bPrime[rowI];

        // Start and end of this row
        fStart = ownStartAddr[rowI];
        fEnd = ownStartAddr[rowI + 1];

        // Accumulate the owner product side
        for (register label curCoeff = fStart; curCoeff < fEnd; curCoeff++)
        {
            // Transpose multiplication.  HJ, 19/Jan/2009
            curX -= lower[curCoeff]*x[upperAddr[curCoeff]];
        }

        // Finish current x
        curX /= diag[rowI];

        // Distribute the neighbour side using current x
        for (register label curCoeff = fStart; curCoeff < fEnd; curCoeff++)
        {
            // Transpose multiplication.  HJ, 19/Jan/2009
            bPrime[upperAddr[curCoeff]] -= upper[curCoeff]*curX;
        }
    }
}


void Foam::coupledGaussSeidelPrecon::reverseSweepTranspose
(
    const lduMatrix& matrix,
    scalarField& x,
    scalarField& bPrime
) const
{
    const scalarField& diag = matrix.diag();
    const scalarField& lower = matrix.lower();
    const scalarField& upper = matrix.upper();

    const labelList& upperAddr = matrix.lduAddr().upperAddr();
    const labelList& ownStartAddr = matrix.lduAddr().ownerStartAddr();

    const label nRows = x.size();
    label fStart, fEnd;

    for (register label rowI = nRows - 1; rowI >= 0; rowI--)
    {
        // lRow is equal to rowI
        scalar& curX = x[rowI];

        // Grab the accumulated neighbour side
        curX = bPrime[rowI];

        // Start and end of this row
        fStart = ownStartAddr[rowI];
        fEnd = ownStartAddr[rowI + 1];

        // Accumulate the owner product side
        for (register label curCoeff = fStart; curCoeff < fEnd; curCoeff++)
        {
            // Transpose multiplication.  HJ, 19/Jan/2009
            curX -= lower[curCoeff]*x[upperAddr[curCoeff]];
        }

        // Finish current x
        curX /= diag[rowI];

        // Distribute the neighbour side using current x
        for (register label curCoeff = fStart; curCoeff < fEnd; curCoeff++)
        {
            // Transpose multiplication.  HJ, 19/Jan/2009
            bPrime[upperAddr[curCoeff]] -= upper[curCoeff]*curX;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coupledGaussSeidelPrecon::coupledGaussSeidelPrecon
(
    const coupledLduMatrix& matrix,
    const PtrList<FieldField<Field, scalar> >& bouCoeffs,
    const PtrList<FieldField<Field, scalar> >& intCoeffs,
    const lduInterfaceFieldPtrsListList& interfaces
)
:
    coupledLduPrecon
    (
        matrix,
        bouCoeffs,
        intCoeffs,
        interfaces
    ),
    mBouCoeffs_(bouCoeffs),
    bPrime_(matrix.size())
{
    // Invert boundary coefficients
    forAll (mBouCoeffs_, rowI)
    {
        mBouCoeffs_[rowI] *= -1;
    }

    // Hook bPrime components
    forAll (matrix_, rowI)
    {
        bPrime_.set(rowI, new scalarField(matrix_[rowI].lduAddr().size(), 0));
    }
}


Foam::coupledGaussSeidelPrecon::coupledGaussSeidelPrecon
(
    const coupledLduMatrix& matrix,
    const PtrList<FieldField<Field, scalar> >& bouCoeffs,
    const PtrList<FieldField<Field, scalar> >& intCoeffs,
    const lduInterfaceFieldPtrsListList& interfaces,
    const dictionary& dict
)
:
    coupledLduPrecon
    (
        matrix,
        bouCoeffs,
        intCoeffs,
        interfaces
    ),
    mBouCoeffs_(bouCoeffs),
    bPrime_(matrix.size())
{
    // Invert boundary coefficients
    forAll (mBouCoeffs_, rowI)
    {
        mBouCoeffs_[rowI] *= -1;
    }

    // Hook bPrime components
    forAll (matrix_, rowI)
    {
        bPrime_.set(rowI, new scalarField(matrix_[rowI].lduAddr().size(), 0));
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::coupledGaussSeidelPrecon::precondition
(
    FieldField<Field, scalar>& x,
    const FieldField<Field, scalar>& b,
    const direction cmpt
) const
{
    // Execute preconditioning
    if (matrix_.diagonal())
    {
        forAll (matrix_, rowI)
        {
            x[rowI] = b[rowI]/matrix_[rowI].diag();
        }
    }
    else if (matrix_.symmetric() || matrix_.asymmetric())
    {
        bPrime_ = b;

        // Parallel boundary update
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
        forAll (matrix_, rowI)
        {
            forwardSweep(matrix_[rowI], x[rowI], bPrime_[rowI]);
        }

        // Reverse sweep
        forAllReverse (matrix_, rowI)
        {
            reverseSweep(matrix_[rowI], x[rowI], bPrime_[rowI]);
        }
    }
}


void Foam::coupledGaussSeidelPrecon::preconditionT
(
    FieldField<Field, scalar>& x,
    const FieldField<Field, scalar>& b,
    const direction cmpt
) const
{
    // Execute preconditioning
    if (matrix_.diagonal())
    {
        forAll (matrix_, rowI)
        {
            x[rowI] = b[rowI]/matrix_[rowI].diag();
        }
    }
    else if (matrix_.symmetric() || matrix_.asymmetric())
    {
        bPrime_ = b;

        // Parallel boundary update
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
        forAll (matrix_, rowI)
        {
            forwardSweepTranspose(matrix_[rowI], x[rowI], bPrime_[rowI]);
        }

        // Reverse sweep
        forAllReverse (matrix_, rowI)
        {
            reverseSweepTranspose(matrix_[rowI], x[rowI], bPrime_[rowI]);
        }
    }
}


// ************************************************************************* //
