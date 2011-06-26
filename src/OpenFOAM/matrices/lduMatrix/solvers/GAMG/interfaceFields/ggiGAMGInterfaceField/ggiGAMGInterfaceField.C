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
    const Pstream::commsTypes commsType
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
    const Pstream::commsTypes commsType
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
            << "Bananas!!!"
            << abort(FatalError);
    }

    forAll(faceCells, elemI)
    {
        result[faceCells[elemI]] -= coeffs[elemI]*zonePnf[elemI];
    }


    // Old treatment
#if(0)
    const labelList& za = ggiInterface_.zoneAddressing();

    forAll(faceCells, elemI)
    {
        result[faceCells[elemI]] -= coeffs[elemI]*zonePnf[za[elemI]];
    }
#endif
}


// ************************************************************************* //
