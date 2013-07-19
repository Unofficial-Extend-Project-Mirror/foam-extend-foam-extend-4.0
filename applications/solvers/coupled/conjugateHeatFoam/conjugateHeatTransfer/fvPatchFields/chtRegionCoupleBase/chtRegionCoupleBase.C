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

Author
    Henrik Rusche, Wikki GmbH.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "chtRcTemperatureFvPatchScalarField.H"
#include "chtRegionCoupleBase.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

chtRegionCoupleBase::chtRegionCoupleBase
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    regionCouplingFvPatchScalarField(p, iF)
{}


chtRegionCoupleBase::chtRegionCoupleBase
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    regionCouplingFvPatchScalarField(p, iF, dict)
{}


chtRegionCoupleBase::chtRegionCoupleBase
(
    const chtRegionCoupleBase& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    regionCouplingFvPatchScalarField(ptf, p, iF, mapper)
{}


chtRegionCoupleBase::chtRegionCoupleBase
(
    const chtRegionCoupleBase& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    regionCouplingFvPatchScalarField(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Return a named shadow patch field
const chtRegionCoupleBase&
chtRegionCoupleBase::shadowPatchField() const
{
    return dynamic_cast<const chtRegionCoupleBase&>
    (
        regionCouplingFvPatchScalarField::shadowPatchField()
    );
}


void chtRegionCoupleBase::initEvaluate
(
    const Pstream::commsTypes commsType
)
{
    updateCoeffs();
}


void
Foam::chtRegionCoupleBase::calcThermalDiffusivity
(
    chtRegionCoupleBase& owner,
    const chtRegionCoupleBase& neighbour
) const
{
    shadowPatchField().calcThermalDiffusivity(owner, neighbour);
}


void
Foam::chtRegionCoupleBase::calcTemperature
(
    chtRcTemperatureFvPatchScalarField& owner,
    const chtRcTemperatureFvPatchScalarField& neighbour,
    const chtRegionCoupleBase& ownerK
) const
{
    shadowPatchField().calcTemperature(owner, neighbour, ownerK);
}


void
Foam::chtRegionCoupleBase::initInterfaceMatrixUpdate
(
    const scalarField&,
    scalarField&,
    const lduMatrix&,
    const scalarField&,
    const direction,
    const Pstream::commsTypes
) const
{
    FatalErrorIn
    (
        "chtRegionCoupleBase::initInterfaceMatrixUpdate"
    )   << abort(FatalError);
}


void
Foam::chtRegionCoupleBase::updateInterfaceMatrix
(
    const scalarField&,
    scalarField&,
    const lduMatrix&,
    const scalarField&,
    const direction,
    const Pstream::commsTypes
) const
{
    FatalErrorIn
    (
        "chtRegionCoupleBase::updateInterfaceMatrix"
    )   << abort(FatalError);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    chtRegionCoupleBase
);

} // End namespace Foam

// ************************************************************************* //
