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
    coupledCholeskyPrecon

Description
    Incomplete Cholesky preconditioning with no fill-in for coupled matrices

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

\*----------------------------------------------------------------------------*/

#include "coupledCholeskyPrecon.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(coupledCholeskyPrecon, 0);

    addToRunTimeSelectionTable
    (
        coupledLduPrecon,
        coupledCholeskyPrecon,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::coupledCholeskyPrecon::calcPreconDiag()
{
    // Precondition the diagonal
    forAll (matrix_, rowI)
    {
        const lduMatrix& rowMatrix = matrix_[rowI];

        preconDiag_.set(rowI, new scalarField(rowMatrix.diag()));
        scalarField& rowPreconDiag = preconDiag_[rowI];

        if (rowMatrix.symmetric())
        {
            const unallocLabelList& upperAddr = rowMatrix.lduAddr().upperAddr();
            const unallocLabelList& lowerAddr = rowMatrix.lduAddr().lowerAddr();

            // Get off-diagonal matrix coefficients
            const scalarField& upper = rowMatrix.upper();

            forAll (upper, coeffI)
            {
                rowPreconDiag[upperAddr[coeffI]] -=
                    sqr(upper[coeffI])/rowPreconDiag[lowerAddr[coeffI]];
            }
        }
        else if (rowMatrix.asymmetric())
        {
            const unallocLabelList& upperAddr = rowMatrix.lduAddr().upperAddr();
            const unallocLabelList& lowerAddr = rowMatrix.lduAddr().lowerAddr();

            // Get off-diagonal matrix coefficients
            const scalarField& upper = rowMatrix.upper();
            const scalarField& lower = rowMatrix.lower();

            forAll (upper, coeffI)
            {
                rowPreconDiag[upperAddr[coeffI]] -=
                    upper[coeffI]*lower[coeffI]/
                    rowPreconDiag[lowerAddr[coeffI]];
            }
        }

        // Invert the diagonal for future use
        forAll (rowPreconDiag, i)
        {
            rowPreconDiag[i] = 1.0/rowPreconDiag[i];
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coupledCholeskyPrecon::coupledCholeskyPrecon
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
    preconDiag_(matrix.size())
{
    calcPreconDiag();
}


Foam::coupledCholeskyPrecon::coupledCholeskyPrecon
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
    preconDiag_(matrix.size())
{
    calcPreconDiag();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::coupledCholeskyPrecon::precondition
(
    FieldField<Field, scalar>& x,
    const FieldField<Field, scalar>& b,
    const direction cmpt
) const
{
    // Cholesky precondition all matrices
    forAll (matrix_, rowI)
    {
        scalarField& rowX = x[rowI];
        const scalarField& rowB = b[rowI];

        const lduMatrix& rowMatrix = matrix_[rowI];
        const scalarField& rowPreconDiag = preconDiag_[rowI];

        forAll(rowX, i)
        {
            rowX[i] = rowB[i]*rowPreconDiag[i];
        }

        if (rowMatrix.symmetric())
        {
            const unallocLabelList& upperAddr = rowMatrix.lduAddr().upperAddr();
            const unallocLabelList& lowerAddr = rowMatrix.lduAddr().lowerAddr();

            // Get off-diagonal matrix coefficients
            const scalarField& upper = rowMatrix.upper();

            forAll (upper, coeffI)
            {
                rowX[upperAddr[coeffI]] -=
                    rowPreconDiag[upperAddr[coeffI]]*
                    upper[coeffI]*rowX[lowerAddr[coeffI]];
            }

            forAllReverse (upper, coeffI)
            {
                rowX[lowerAddr[coeffI]] -=
                    rowPreconDiag[lowerAddr[coeffI]]*
                    upper[coeffI]*rowX[upperAddr[coeffI]];
            }
        }
        else if (rowMatrix.asymmetric())
        {
            const unallocLabelList& upperAddr = rowMatrix.lduAddr().upperAddr();
            const unallocLabelList& lowerAddr = rowMatrix.lduAddr().lowerAddr();
            const unallocLabelList& losortAddr =
                rowMatrix.lduAddr().losortAddr();

            // Get off-diagonal matrix coefficients
            const scalarField& upper = rowMatrix.upper();
            const scalarField& lower = rowMatrix.lower();

            label losortCoeff;

            forAll (lower, coeffI)
            {
                losortCoeff = losortAddr[coeffI];

                rowX[upperAddr[losortCoeff]] -=
                    rowPreconDiag[upperAddr[losortCoeff]]*
                    lower[losortCoeff]*rowX[lowerAddr[losortCoeff]];
            }

            forAllReverse (upper, coeffI)
            {
                rowX[lowerAddr[coeffI]] -=
                    rowPreconDiag[lowerAddr[coeffI]]*
                    upper[coeffI]*rowX[upperAddr[coeffI]];
            }
        }
    }
}


void Foam::coupledCholeskyPrecon::preconditionT
(
    FieldField<Field, scalar>& x,
    const FieldField<Field, scalar>& b,
    const direction cmpt
) const
{
    // Cholesky precondition all matrices
    forAll (matrix_, rowI)
    {
        scalarField& rowX = x[rowI];
        const scalarField& rowB = b[rowI];

        const lduMatrix& rowMatrix = matrix_[rowI];
        const scalarField& rowPreconDiag = preconDiag_[rowI];

        forAll(rowX, i)
        {
            rowX[i] = rowB[i]*rowPreconDiag[i];
        }

        if (rowMatrix.symmetric())
        {
            const unallocLabelList& upperAddr = rowMatrix.lduAddr().upperAddr();
            const unallocLabelList& lowerAddr = rowMatrix.lduAddr().lowerAddr();

            // Get off-diagonal matrix coefficients
            const scalarField& upper = rowMatrix.upper();

            forAll (upper, coeffI)
            {
                // For symmetric matrix, there is no change.  HJ, 19/Jan/2009
                rowX[upperAddr[coeffI]] -=
                    rowPreconDiag[upperAddr[coeffI]]*
                    upper[coeffI]*rowX[lowerAddr[coeffI]];
            }

            forAllReverse (upper, coeffI)
            {
                // For symmetric matrix, there is no change.  HJ, 19/Jan/2009
                rowX[lowerAddr[coeffI]] -=
                    rowPreconDiag[lowerAddr[coeffI]]*
                    upper[coeffI]*rowX[upperAddr[coeffI]];
            }
        }
        else if (rowMatrix.asymmetric())
        {
            const unallocLabelList& upperAddr = rowMatrix.lduAddr().upperAddr();
            const unallocLabelList& lowerAddr = rowMatrix.lduAddr().lowerAddr();
            const unallocLabelList& losortAddr =
                rowMatrix.lduAddr().losortAddr();

            // Get off-diagonal matrix coefficients
            const scalarField& upper = rowMatrix.upper();
            const scalarField& lower = rowMatrix.lower();

            label losortCoeff;

            forAll (lower, coeffI)
            {
                // Transpose multiplication.  HJ, 19/Jan/2009
                rowX[upperAddr[coeffI]] -=
                    rowPreconDiag[upperAddr[coeffI]]*
                    upper[coeffI]*rowX[lowerAddr[coeffI]];
            }

            forAllReverse (upper, coeffI)
            {
                losortCoeff = losortAddr[coeffI];

                // Transpose multiplication.  HJ, 19/Jan/2009
                rowX[lowerAddr[losortCoeff]] -=
                    rowPreconDiag[lowerAddr[losortCoeff]]*
                    lower[coeffI]*rowX[upperAddr[losortCoeff]];
            }
        }
    }
}


// ************************************************************************* //
