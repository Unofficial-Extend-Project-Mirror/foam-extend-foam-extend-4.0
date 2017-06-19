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
    CholeskyPrecon

Description
    Incmplete Cholesky preconditioning with no fill-in

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "CholeskyPrecon.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(CholeskyPrecon, 0);

    lduPreconditioner::
        addsymMatrixConstructorToTable<CholeskyPrecon>
        addCholeskyPreconditionerSymMatrixConstructorToTable_;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::CholeskyPrecon::calcPreconDiag()
{
    // Precondition the diagonal
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

        forAll (upper, coeffI)
        {
            preconDiag_[upperAddr[coeffI]] -=
                sqr(upper[coeffI])/preconDiag_[lowerAddr[coeffI]];
        }
    }

    // Invert the diagonal for future use
    forAll (preconDiag_, i)
    {
        preconDiag_[i] = 1.0/preconDiag_[i];
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::CholeskyPrecon::CholeskyPrecon
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


Foam::CholeskyPrecon::CholeskyPrecon
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

Foam::CholeskyPrecon::~CholeskyPrecon()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::CholeskyPrecon::precondition
(
    scalarField& x,
    const scalarField& b,
    const direction cmpt
) const
{
    if (matrix_.asymmetric())
    {
        FatalErrorIn
        (
            "void CholeskyPrecon::precondition\n"
            "(\n"
            "    scalarField& x,\n"
            "    const scalarField& b,\n"
            "    const direction cmpt\n"
            ") const"
        )   << "Calling CholeskyPrecon on an assymetric matrix.  "
            << "Please use ILU0 instead"
            << abort(FatalError);
    }

    // In order to properly execute parallel preconditioning, re-use
    // x to zero and execute coupled boundary update
    // HJ, 19/Jun/2017

    x = 0;

    // Coupled boundary update
    {
        matrix_.initMatrixInterfaces
        (
            coupleBouCoeffs_,
            interfaces_,
            x,
            x,               // put result into x
            cmpt,
            false
        );

        matrix_.updateMatrixInterfaces
        (
            coupleBouCoeffs_,
            interfaces_,
            x,
            x,               // put result into x
            cmpt,
            false
        );
    }

    // Multiply with inverse diag to precondition
    x *= preconDiag_;

    // Diagonal block: note +=
    {
        scalar* __restrict__ xPtr = x.begin();

        const scalar* __restrict__ preconDiagPtr = preconDiag_.begin();

        const scalar* __restrict__ bPtr = b.begin();

        const label nRows = x.size();

        for (register label rowI = 0; rowI < nRows; rowI++)
        {
            xPtr[rowI] += bPtr[rowI]*preconDiagPtr[rowI];
        }
    }

    if (matrix_.symmetric())
    {
        scalar* __restrict__ xPtr = x.begin();

        // Addressing
        const label* const __restrict__ uPtr =
            matrix_.lduAddr().upperAddr().begin();

        const label* const __restrict__ lPtr =
            matrix_.lduAddr().lowerAddr().begin();

        // Coeffs
        const scalar* __restrict__ preconDiagPtr = preconDiag_.begin();

        const scalar* const __restrict__ upperPtr = matrix_.upper().begin();

        const label nCoeffs = matrix_.upper().size();

        // Forward sweep
        for (register label coeffI = 0; coeffI < nCoeffs; coeffI++)
        {
            xPtr[uPtr[coeffI]] -=
                preconDiagPtr[uPtr[coeffI]]*
                upperPtr[coeffI]*xPtr[lPtr[coeffI]];
        }

        // Reverse sweep
        for (register label coeffI = nCoeffs - 1; coeffI >= 0; coeffI--)
        {
            xPtr[lPtr[coeffI]] -=
                preconDiagPtr[lPtr[coeffI]]*
                upperPtr[coeffI]*xPtr[uPtr[coeffI]];
        }
    }
}


// ************************************************************************* //
