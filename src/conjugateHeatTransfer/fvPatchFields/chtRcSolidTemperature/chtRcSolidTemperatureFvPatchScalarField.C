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

#include "chtRcSolidTemperatureFvPatchScalarField.H"
#include "chtRegionCoupleBase.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "fvMatrices.H"
#include "basicThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::chtRcSolidTemperatureFvPatchScalarField::
chtRcSolidTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    chtRcTemperatureFvPatchScalarField(p, iF)
{}


Foam::chtRcSolidTemperatureFvPatchScalarField::
chtRcSolidTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    chtRcTemperatureFvPatchScalarField(p, iF, dict)
{}


Foam::chtRcSolidTemperatureFvPatchScalarField::
chtRcSolidTemperatureFvPatchScalarField
(
    const chtRcSolidTemperatureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    chtRcTemperatureFvPatchScalarField(ptf, p, iF, mapper)
{}


Foam::chtRcSolidTemperatureFvPatchScalarField::
chtRcSolidTemperatureFvPatchScalarField
(
    const chtRcSolidTemperatureFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    chtRcTemperatureFvPatchScalarField(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// Return neighbour field given internal cell data
Foam::tmp<Foam::scalarField>
Foam::chtRcSolidTemperatureFvPatchScalarField::patchNeighbourField() const
{
    return regionCouplingFvPatchScalarField::patchNeighbourField("T");
}


// Return a named shadow patch field
const Foam::chtRcTemperatureFvPatchScalarField&
Foam::chtRcSolidTemperatureFvPatchScalarField::shadowPatchField() const
{
    return refCast<const chtRcTemperatureFvPatchScalarField>
    (
        regionCouplingFvPatchScalarField::shadowPatchField()
    );
}


void Foam::chtRcSolidTemperatureFvPatchScalarField::evaluate
(
    const Pstream::commsTypes
)
{
    fvPatchScalarField::evaluate();
}


void Foam::chtRcSolidTemperatureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const chtRegionCoupleBase& K =
        refCast<const chtRegionCoupleBase>
        (
            lookupPatchField<volScalarField, scalar>(kName())
        );

    *this == K.calcTemperature(*this, shadowPatchField(), K);

    fvPatchScalarField::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

makePatchTypeField
(
    fvPatchScalarField,
    chtRcSolidTemperatureFvPatchScalarField
);

} // End namespace Foam

// ************************************************************************* //
