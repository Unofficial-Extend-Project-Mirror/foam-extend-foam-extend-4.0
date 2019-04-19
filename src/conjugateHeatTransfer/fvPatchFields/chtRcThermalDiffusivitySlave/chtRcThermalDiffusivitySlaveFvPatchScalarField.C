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

Author
    Henrik Rusche, Wikki GmbH.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "chtRcThermalDiffusivitySlaveFvPatchScalarField.H"
#include "chtRcTemperatureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "basicThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::chtRcThermalDiffusivitySlaveFvPatchScalarField::
chtRcThermalDiffusivitySlaveFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    chtRegionCoupleBase(p, iF)
{}


Foam::chtRcThermalDiffusivitySlaveFvPatchScalarField::
chtRcThermalDiffusivitySlaveFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    chtRegionCoupleBase(p, iF, dict)
{}



Foam::chtRcThermalDiffusivitySlaveFvPatchScalarField::
chtRcThermalDiffusivitySlaveFvPatchScalarField
(
    const chtRcThermalDiffusivitySlaveFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    chtRegionCoupleBase(ptf, p, iF, mapper)
{}


Foam::chtRcThermalDiffusivitySlaveFvPatchScalarField::
chtRcThermalDiffusivitySlaveFvPatchScalarField
(
    const chtRcThermalDiffusivitySlaveFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    chtRegionCoupleBase(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::chtRegionCoupleBase&
Foam::chtRcThermalDiffusivitySlaveFvPatchScalarField::shadowPatchField() const
{
    return dynamic_cast<const chtRegionCoupleBase&>
    (
        regionCouplingFvPatchScalarField::shadowPatchField()
    );
}


void Foam::chtRcThermalDiffusivitySlaveFvPatchScalarField::evaluate
(
    const Pstream::commsTypes
)
{
    fvPatchScalarField::evaluate();
}


void Foam::chtRcThermalDiffusivitySlaveFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    scalarField& k = *this;

    if
    (
        dimensionedInternalField().dimensions()
     == dimensionSet(1, -1, -1, 0, 0, 0, 0)
    )
    {
        const label patchi = patch().index();

        const basicThermo& thermo = db().lookupObject<basicThermo>
        (
            "thermophysicalProperties"
        );

        const chtRcTemperatureFvPatchScalarField& h =
            refCast<const chtRcTemperatureFvPatchScalarField>
            (
                thermo.h().boundaryField()[patchi]
            );

        const scalarField Tw =
            lookupPatchField<volScalarField, scalar>("T");

        k = calcThermalDiffusivity(*this, shadowPatchField(), h)
            /thermo.Cp(Tw, patchi);
    }
    else
    {
        const chtRcTemperatureFvPatchScalarField& Tw =
            refCast<const chtRcTemperatureFvPatchScalarField>
            (
                lookupPatchField<volScalarField, scalar>("T")
            );

        k = calcThermalDiffusivity(*this, shadowPatchField(), Tw);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    makePatchTypeField
    (
        fvPatchScalarField,
        chtRcThermalDiffusivitySlaveFvPatchScalarField
    );

} // End namespace Foam


// ************************************************************************* //
