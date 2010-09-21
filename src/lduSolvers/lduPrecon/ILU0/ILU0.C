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
    if (matrix_.asymmetric())
    {
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
    const direction
) const
{
    forAll(x, i)
    {
        x[i] = b[i]*preconDiag_[i];
    }

    if (matrix_.asymmetric())
    {
        const unallocLabelList& upperAddr = matrix_.lduAddr().upperAddr();
        const unallocLabelList& lowerAddr = matrix_.lduAddr().lowerAddr();
        const unallocLabelList& losortAddr = matrix_.lduAddr().losortAddr();

        // Get off-diagonal matrix coefficients
        const scalarField& upper = matrix_.upper();
        const scalarField& lower = matrix_.lower();

        label losortCoeff;

        forAll (lower, coeffI)
        {
            losortCoeff = losortAddr[coeffI];

            x[upperAddr[losortCoeff]] -=
                preconDiag_[upperAddr[losortCoeff]]*
                lower[losortCoeff]*x[lowerAddr[losortCoeff]];
        }

        forAllReverse (upper, coeffI)
        {
            x[lowerAddr[coeffI]] -=
                preconDiag_[lowerAddr[coeffI]]*
                upper[coeffI]*x[upperAddr[coeffI]];
        }
    }
}


void Foam::ILU0::preconditionT
(
    scalarField& x,
    const scalarField& b,
    const direction cmpt
) const
{
    forAll(x, i)
    {
        x[i] = b[i]*preconDiag_[i];
    }

    if (matrix_.asymmetric())
    {
        const unallocLabelList& upperAddr = matrix_.lduAddr().upperAddr();
        const unallocLabelList& lowerAddr = matrix_.lduAddr().lowerAddr();
        const unallocLabelList& losortAddr = matrix_.lduAddr().losortAddr();

        // Get off-diagonal matrix coefficients
        const scalarField& upper = matrix_.upper();
        const scalarField& lower = matrix_.lower();

        label losortCoeff;

        forAll (lower, coeffI)
        {
            // Transpose multiplication.  HJ, 19/Jan/2009
            x[upperAddr[coeffI]] -=
                preconDiag_[upperAddr[coeffI]]*
                upper[coeffI]*x[lowerAddr[coeffI]];
        }

        forAllReverse (upper, coeffI)
        {
            losortCoeff = losortAddr[coeffI];

            // Transpose multiplication.  HJ, 19/Jan/2009
            x[lowerAddr[losortCoeff]] -=
                preconDiag_[lowerAddr[losortCoeff]]*
                lower[losortCoeff]*x[upperAddr[losortCoeff]];
        }
    }
}


// ************************************************************************* //
