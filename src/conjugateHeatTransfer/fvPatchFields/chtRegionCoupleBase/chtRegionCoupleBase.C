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

#include "chtRcTemperatureFvPatchScalarField.H"
#include "chtRegionCoupleBase.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "fvMatrices.H"
#include "magLongDelta.H"
#include "basicThermo.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::chtRegionCoupleBase::chtRegionCoupleBase
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    regionCouplingFvPatchScalarField(p, iF)
{}


Foam::chtRegionCoupleBase::chtRegionCoupleBase
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    regionCouplingFvPatchScalarField(p, iF, dict)
{}


Foam::chtRegionCoupleBase::chtRegionCoupleBase
(
    const chtRegionCoupleBase& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    regionCouplingFvPatchScalarField(ptf, p, iF, mapper)
{}


Foam::chtRegionCoupleBase::chtRegionCoupleBase
(
    const chtRegionCoupleBase& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    regionCouplingFvPatchScalarField(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Return a named shadow patch field
const Foam::chtRegionCoupleBase&
Foam::chtRegionCoupleBase::shadowPatchField() const
{
    return dynamic_cast<const chtRegionCoupleBase&>
    (
        regionCouplingFvPatchScalarField::shadowPatchField()
    );
}


// Return a named shadow patch field
Foam::tmp<Foam::scalarField> Foam::chtRegionCoupleBase::forig() const
{
    if
    (
        dimensionedInternalField().dimensions()
     == dimensionSet(1, -1, -1, 0, 0)
    )
    {
        const fvPatch& p = patch();

        const basicThermo& thermo = db().lookupObject<basicThermo>
        (
            "thermophysicalProperties"
        );

        const scalarField Tw =
            p.lookupPatchField<volScalarField, scalar>("T");

        return originalPatchField()*thermo.Cp(Tw, p.index());
    }
    else
    {
        return originalPatchField();
    }
}


// Return a named shadow patch field
Foam::tmp<Foam::scalarField> Foam::chtRegionCoupleBase::korig() const
{
    const fvPatch& p = patch();
    const magLongDelta& mld = magLongDelta::New(p.boundaryMesh().mesh());

    return forig()/(1 - p.weights())/mld.magDelta(p.index());
}


// Return a named shadow patch field
Foam::tmp<Foam::scalarField> Foam::chtRegionCoupleBase::kw() const
{
    if
    (
        dimensionedInternalField().dimensions()
     == dimensionSet(1, -1, -1, 0, 0, 0, 0)
    )
    {
        const fvPatch& p = patch();

        const basicThermo& thermo = db().lookupObject<basicThermo>
        (
            "thermophysicalProperties"
        );

        const scalarField Tw =
            p.lookupPatchField<volScalarField, scalar>("T");

        return *this*thermo.Cp(Tw, p.index());
    }
    else
    {
        return *this;
    }
}


Foam::tmp<Foam::scalarField>
Foam::chtRegionCoupleBase::calcThermalDiffusivity
(
    const chtRegionCoupleBase& owner,
    const chtRegionCoupleBase& neighbour,
    const chtRcTemperatureFvPatchScalarField& TwOwn
) const
{
    return shadowPatchField().calcThermalDiffusivity(owner, neighbour, TwOwn);
}


Foam::tmp<Foam::scalarField>
Foam::chtRegionCoupleBase::calcTemperature
(
    const chtRcTemperatureFvPatchScalarField& TwOwn,
    const chtRcTemperatureFvPatchScalarField& neighbour,
    const chtRegionCoupleBase& ownerK
) const
{
    return shadowPatchField().calcTemperature(TwOwn, neighbour, ownerK);
}


void Foam::chtRegionCoupleBase::initEvaluate
(
    const Pstream::commsTypes commsType
)
{
    updateCoeffs();
}


void
Foam::chtRegionCoupleBase::initInterfaceMatrixUpdate
(
    const scalarField&,
    scalarField&,
    const lduMatrix&,
    const scalarField&,
    const direction,
    const Pstream::commsTypes,
    const bool switchToLhs
) const
{
    FatalErrorIn
    (
        "chtRegionCoupleBase::initInterfaceMatrixUpdate"
    )   << "Undefined function: this patch field cannot be used "
        << "on active variables"
        << abort(FatalError);
}


void
Foam::chtRegionCoupleBase::updateInterfaceMatrix
(
    const scalarField&,
    scalarField&,
    const lduMatrix&,
    const scalarField&,
    const direction,
    const Pstream::commsTypes,
    const bool switchToLhs
) const
{
    FatalErrorIn
    (
        "chtRegionCoupleBase::updateInterfaceMatrix"
    )   << "Undefined function: this patch field cannot be used "
        << "on active variables"
        << abort(FatalError);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

makePatchTypeField
(
    fvPatchScalarField,
    chtRegionCoupleBase
);

} // End namespace Foam

// ************************************************************************* //
