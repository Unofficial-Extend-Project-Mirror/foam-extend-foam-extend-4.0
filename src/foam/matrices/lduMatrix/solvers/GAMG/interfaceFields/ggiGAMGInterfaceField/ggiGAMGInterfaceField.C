/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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

#include "ggiGAMGInterfaceField.H"
#include "addToRunTimeSelectionTable.H"
#include "lduMatrix.H"

#include "OSspecific.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ggiGAMGInterfaceField, 0);
    addToRunTimeSelectionTable
    (
        GAMGInterfaceField,
        ggiGAMGInterfaceField,
        lduInterface
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ggiGAMGInterfaceField::ggiGAMGInterfaceField
(
    const GAMGInterface& GAMGCp,
    const lduInterfaceField& fineInterface
)
:
    GAMGInterfaceField(GAMGCp, fineInterface),
    ggiInterface_(refCast<const ggiGAMGInterface>(GAMGCp)),
    doTransform_(false),
    rank_(0),
    transferBuffer_()
{
    const ggiLduInterfaceField& p =
        refCast<const ggiLduInterfaceField>(fineInterface);

    doTransform_ = p.doTransform();
    rank_ = p.rank();
}


// * * * * * * * * * * * * * * * * Desstructor * * * * * * * * * * * * * * * //

Foam::ggiGAMGInterfaceField::~ggiGAMGInterfaceField()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::ggiGAMGInterfaceField::initInterfaceMatrixUpdate
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


void Foam::ggiGAMGInterfaceField::updateInterfaceMatrix
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
    scalarField zonePnf =
        ggiInterface_.internalFieldTransfer(commsType, psiInternal);
    transformCoupleField(zonePnf, cmpt);

    const unallocLabelList& faceCells = ggiInterface_.faceCells();

    // New treatment.  HJ, 26/Jun/2011
    if (zonePnf.size() != faceCells.size())
    {
        FatalErrorIn("ggiGAMGInterfaceField::updateInterfaceMatrix")
            << "Error in interface update: incorrect size of zone fields" << nl
            << "Field size = " << zonePnf.size()
            << " Zone size = " << faceCells.size()
            << abort(FatalError);
    }

    if (switchToLhs)
    {
        forAll(faceCells, elemI)
        {
            result[faceCells[elemI]] += coeffs[elemI]*zonePnf[elemI];
        }
    }
    else
    {
        forAll(faceCells, elemI)
        {
            result[faceCells[elemI]] -= coeffs[elemI]*zonePnf[elemI];
        }
    }
}


// ************************************************************************* //
