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

#include "chtRcThermalDiffusivitySlaveFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

chtRcThermalDiffusivitySlaveFvPatchScalarField::chtRcThermalDiffusivitySlaveFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    chtRegionCoupleBase(p, iF)
{}


chtRcThermalDiffusivitySlaveFvPatchScalarField::chtRcThermalDiffusivitySlaveFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    chtRegionCoupleBase(p, iF, dict)
{}


chtRcThermalDiffusivitySlaveFvPatchScalarField::chtRcThermalDiffusivitySlaveFvPatchScalarField
(
    const chtRcThermalDiffusivitySlaveFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    chtRegionCoupleBase(ptf, p, iF, mapper)
{}


chtRcThermalDiffusivitySlaveFvPatchScalarField::chtRcThermalDiffusivitySlaveFvPatchScalarField
(
    const chtRcThermalDiffusivitySlaveFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    chtRegionCoupleBase(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const chtRegionCoupleBase&
chtRcThermalDiffusivitySlaveFvPatchScalarField::shadowPatchField() const
{
    return dynamic_cast<const chtRegionCoupleBase&>
    (
        regionCouplingFvPatchScalarField::shadowPatchField()
    );
}


void chtRcThermalDiffusivitySlaveFvPatchScalarField::evaluate
(
    const Pstream::commsTypes
)
{
    fvPatchScalarField::evaluate();
}


void chtRcThermalDiffusivitySlaveFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    shadowPatchField().calcThermalDiffusivity(*this, shadowPatchField());
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    chtRcThermalDiffusivitySlaveFvPatchScalarField
);

} // End namespace Foam

// ************************************************************************* //
