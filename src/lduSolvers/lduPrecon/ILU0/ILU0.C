/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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
    ILU0

Description
    Incmplete Cholesky preconditioning with no fill-in

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "ILU0.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ILU0, 0);

    lduPreconditioner::
        addasymMatrixConstructorToTable<ILU0>
        addILU0ditionerAsymMatrixConstructorToTable_;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::ILU0::calcPreconDiag()
{
    // Precondition the diagonal

    // Do coupled interfaces

    // Note: ordering of ILU coefficients on the coupled boundary is
    // out of sequence.  Formally, this can be fixed by visiting the coupled
    // boundary from the cell loop, but HJ thinks this does not
    // make a difference
    // HJ and VV, 19/Jun/2017
    if (!matrix_.diagonal())
    {
        // Do coupled interfaces
        forAll (interfaces_, patchI)
        {
            if (interfaces_.set(patchI))
            {
                // Get face-cells addressing
                const unallocLabelList& fc =
                    interfaces_[patchI].coupledInterface().faceCells();

                // Get interface coefficiens
                const scalarField& bouCoeffs = coupleBouCoeffs_[patchI];
                const scalarField& intCoeffs = coupleIntCoeffs_[patchI];

                forAll (fc, coeffI)
                {
                    // Note: change of the sign compared to main loop below
                    // This is because lower = -intCoeffs
                    // HJ and VV, 19/Jun/2017
                    // Note: sign fixed by HJ, 19/Jun/2017
                    preconDiag_[fc[coeffI]] +=
                        bouCoeffs[coeffI]*intCoeffs[coeffI]/
                        preconDiag_[fc[coeffI]];
                }
            }
        }

        // Do core matrix

        const unallocLabelList& upperAddr = matrix_.lduAddr().upperAddr();
        const unallocLabelList& lowerAddr = matrix_.lduAddr().lowerAddr();

        // Get off-diagonal matrix coefficients
        const scalarField& upper = matrix_.upper();
        const scalarField& lower = matrix_.lower();

        forAll (upper, coeffI)
        {
            preconDiag_[upperAddr[coeffI]] -=
                upper[coeffI]*lower[coeffI]/preconDiag_[lowerAddr[coeffI]];
        }
    }

    // Invert the diagonal for future use
    forAll (preconDiag_, i)
    {
        preconDiag_[i] = 1.0/preconDiag_[i];
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ILU0::ILU0
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
    preconDiag_(matrix_.diag())
{
    calcPreconDiag();
}


Foam::ILU0::ILU0
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
    preconDiag_(matrix_.diag())
{
    calcPreconDiag();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ILU0::~ILU0()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::ILU0::precondition
(
    scalarField& x,
    const scalarField& b,
    const direction cmpt
) const
{
    if (matrix_.symmetric())
    {
        FatalErrorIn
        (
            "void ILU0::precondition\n"
            "(\n"
            "    scalarField& x,\n"
            "    const scalarField& b,\n"
            "    const direction cmpt\n"
            ") const"
        )   << "Calling ILU0 on a symetric matrix.  "
            << "Please use CholeskyPrecon instead"
            << abort(FatalError);
    }

    // Note: coupled boundary updated is not needed because x is zero
    // HJ and VV, 19/Jun/2017

    // Diagonal block
    {
        scalar* __restrict__ xPtr = x.begin();

        const scalar* __restrict__ preconDiagPtr = preconDiag_.begin();

        const scalar* __restrict__ bPtr = b.begin();

        const label nRows = x.size();

        // Note: multiplication over-write x: no need to initialise
        // HJ, and VV, 19/Jun/2017
        for (label rowI = 0; rowI < nRows; rowI++)
        {
            xPtr[rowI] = bPtr[rowI]*preconDiagPtr[rowI];
        }
    }

    if (matrix_.asymmetric())
    {
        const unallocLabelList& upperAddr = matrix_.lduAddr().upperAddr();
        const unallocLabelList& lowerAddr = matrix_.lduAddr().lowerAddr();
        const unallocLabelList& losortAddr = matrix_.lduAddr().losortAddr();

        // Get off-diagonal matrix coefficients
        const scalarField& upper = matrix_.upper();
        const scalarField& lower = matrix_.lower();

        label losortIndex;

        forAll (lower, coeffI)
        {
            losortIndex = losortAddr[coeffI];

            x[upperAddr[losortIndex]] -=
                preconDiag_[upperAddr[losortIndex]]*
                lower[losortIndex]*x[lowerAddr[losortIndex]];
        }

        forAllReverse (upper, coeffI)
        {
            x[lowerAddr[coeffI]] -=
                preconDiag_[lowerAddr[coeffI]]*
                upper[coeffI]*x[upperAddr[coeffI]];
        }
/*
        // Parallel preconditioning
        // HJ, 19/Jun/2017

        scalarField xCorr(x.size(), 0);

        // Coupled boundary update
        {
            matrix_.initMatrixInterfaces
            (
                coupleBouCoeffs_,
                interfaces_,
                x,
                xCorr,               // put result into xCorr
                cmpt,
                false
            );

            matrix_.updateMatrixInterfaces
            (
                coupleBouCoeffs_,
                interfaces_,
                x,
                xCorr,               // put result into xCorr
                cmpt,
                false
            );

            // Multiply with inverse diag to precondition
            x += xCorr*preconDiag_;
        }
*/
    }
}


void Foam::ILU0::preconditionT
(
    scalarField& x,
    const scalarField& b,
    const direction cmpt
) const
{
    if (matrix_.symmetric())
    {
        FatalErrorIn
        (
            "void ILU0::precondition\n"
            "(\n"
            "    scalarField& x,\n"
            "    const scalarField& b,\n"
            "    const direction cmpt\n"
            ") const"
        )   << "Calling ILU0 on a symetric matrix.  "
            << "Please use CholeskyPrecon instead"
            << abort(FatalError);
    }

    // Note: coupled boundary updated is not needed because x is zero
    // HJ and VV, 19/Jun/2017

    // Diagonal block
    {
        scalar* __restrict__ xPtr = x.begin();

        const scalar* __restrict__ preconDiagPtr = preconDiag_.begin();

        const scalar* __restrict__ bPtr = b.begin();

        const label nRows = x.size();

        // Note: multiplication over-write x: no need to initialise
        // HJ, and VV, 19/Jun/2017
        for (label rowI = 0; rowI < nRows; rowI++)
        {
            xPtr[rowI] = bPtr[rowI]*preconDiagPtr[rowI];
        }
    }

    if (matrix_.asymmetric())
    {
        const unallocLabelList& upperAddr = matrix_.lduAddr().upperAddr();
        const unallocLabelList& lowerAddr = matrix_.lduAddr().lowerAddr();
        const unallocLabelList& losortAddr = matrix_.lduAddr().losortAddr();

        // Get off-diagonal matrix coefficients
        const scalarField& upper = matrix_.upper();
        const scalarField& lower = matrix_.lower();

        label losortIndex;

        forAll (lower, coeffI)
        {
            // Transpose multiplication.  HJ, 19/Jan/2009
            x[upperAddr[coeffI]] -=
                preconDiag_[upperAddr[coeffI]]*
                upper[coeffI]*x[lowerAddr[coeffI]];
        }

        forAllReverse (upper, coeffI)
        {
            losortIndex = losortAddr[coeffI];

            // Transpose multiplication.  HJ, 19/Jan/2009
            x[lowerAddr[losortIndex]] -=
                preconDiag_[lowerAddr[losortIndex]]*
                lower[losortIndex]*x[upperAddr[losortIndex]];
        }
    }
}


// ************************************************************************* //
