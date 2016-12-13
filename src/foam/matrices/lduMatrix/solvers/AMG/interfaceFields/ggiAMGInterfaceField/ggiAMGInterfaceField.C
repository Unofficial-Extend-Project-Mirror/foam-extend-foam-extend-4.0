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

\*---------------------------------------------------------------------------*/

#include "ggiAMGInterfaceField.H"
#include "addToRunTimeSelectionTable.H"
#include "lduMatrix.H"

#include "OSspecific.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ggiAMGInterfaceField, 0);
    addToRunTimeSelectionTable
    (
        AMGInterfaceField,
        ggiAMGInterfaceField,
        lduInterface
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ggiAMGInterfaceField::ggiAMGInterfaceField
(
    const AMGInterface& AMGCp,
    const lduInterfaceField& fineInterface
)
:
    AMGInterfaceField(AMGCp, fineInterface),
    ggiInterface_(refCast<const ggiAMGInterface>(AMGCp)),
    doTransform_(false),
    rank_(0)
{
    const ggiLduInterfaceField& p =
        refCast<const ggiLduInterfaceField>(fineInterface);

    doTransform_ = p.doTransform();
    rank_ = p.rank();
}


// * * * * * * * * * * * * * * * * Desstructor * * * * * * * * * * * * * * * //

Foam::ggiAMGInterfaceField::~ggiAMGInterfaceField()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::ggiAMGInterfaceField::initInterfaceMatrixUpdate
(
    const scalarField& psiInternal,
    scalarField& result,
    const lduMatrix&,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes commsType,
    const bool switchToLhs
) const
{
    // This must have a reduce in it.  HJ, 15/May/2009
    ggiInterface_.initInternalFieldTransfer(commsType, psiInternal);
}


void Foam::ggiAMGInterfaceField::updateInterfaceMatrix
(
    const scalarField& psiInternal,
    scalarField& result,
    const lduMatrix&,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes commsType,
    const bool switchToLhs
) const
{
    // Get expanded data to zone size.  No global reduce allowed
    // HJ, 15/May/2009
    scalarField pnf =
        ggiInterface_.internalFieldTransfer(commsType, psiInternal);
    transformCoupleField(pnf, cmpt);

    const unallocLabelList& faceCells = ggiInterface_.faceCells();

    // New treatment.  HJ, 26/Jun/2011
    if (pnf.size() != faceCells.size())
    {
        FatalErrorIn("ggiAMGInterfaceField::updateInterfaceMatrix")
            << "Error in interface update: incorrect size of zone fields" << nl
            << "Field size = " << pnf.size()
            << " Zone size = " << faceCells.size()
            << abort(FatalError);
    }

    if (switchToLhs)
    {
        forAll(faceCells, elemI)
        {
            result[faceCells[elemI]] += coeffs[elemI]*pnf[elemI];
        }
    }
    else
    {
        forAll(faceCells, elemI)
        {
            result[faceCells[elemI]] -= coeffs[elemI]*pnf[elemI];
        }
    }
}


// ************************************************************************* //
