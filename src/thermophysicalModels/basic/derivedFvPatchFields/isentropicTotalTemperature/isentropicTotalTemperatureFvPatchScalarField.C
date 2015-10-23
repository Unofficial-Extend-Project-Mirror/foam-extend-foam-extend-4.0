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

#include "isentropicTotalTemperatureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "basicThermo.H"
#include "isentropicTotalPressureFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::isentropicTotalTemperatureFvPatchScalarField::
isentropicTotalTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    pName_("p"),
    T0_(p.size(), 0.0)
{}


Foam::isentropicTotalTemperatureFvPatchScalarField::
isentropicTotalTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    pName_(dict.lookupOrDefault<word>("p", "p")),
    T0_("T0", dict, p.size())
{
    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        fvPatchField<scalar>::operator=(T0_);
    }
}


Foam::isentropicTotalTemperatureFvPatchScalarField::
isentropicTotalTemperatureFvPatchScalarField
(
    const isentropicTotalTemperatureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    pName_(ptf.pName_),
    T0_(ptf.T0_, mapper)
{}


Foam::isentropicTotalTemperatureFvPatchScalarField::
isentropicTotalTemperatureFvPatchScalarField
(
    const isentropicTotalTemperatureFvPatchScalarField& ptf
)
:
    fixedValueFvPatchScalarField(ptf),
    pName_(ptf.pName_),
    T0_(ptf.T0_)
{}


Foam::isentropicTotalTemperatureFvPatchScalarField::
isentropicTotalTemperatureFvPatchScalarField
(
    const isentropicTotalTemperatureFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    pName_(ptf.pName_),
    T0_(ptf.T0_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::isentropicTotalTemperatureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchScalarField::autoMap(m);
    T0_.autoMap(m);
}


void Foam::isentropicTotalTemperatureFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchScalarField::rmap(ptf, addr);

    const isentropicTotalTemperatureFvPatchScalarField& tiptf =
        refCast<const isentropicTotalTemperatureFvPatchScalarField>(ptf);

    T0_.rmap(tiptf.T0_, addr);
}


void Foam::isentropicTotalTemperatureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Get pressure and temperature
    const scalarField& T = *this;

    const fvPatchScalarField& pp =
        patch().lookupPatchField<volScalarField, scalar>(pName_);

    const isentropicTotalPressureFvPatchScalarField& p =
        refCast<const isentropicTotalPressureFvPatchScalarField>(pp);

    const basicThermo& thermo =
        db().lookupObject<basicThermo>("thermophysicalProperties");

    scalarField gamma =
        thermo.Cp(T, patch().index())/thermo.Cv(T, patch().index());

    operator==(T0_*pow(p/p.p0(), (gamma - 1)/gamma));

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::isentropicTotalTemperatureFvPatchScalarField::updateCoeffs
(
    const vectorField& Up
)
{
    updateCoeffs();
}


void Foam::isentropicTotalTemperatureFvPatchScalarField::write
(
    Ostream& os
) const
{
    fixedValueFvPatchScalarField::write(os);
    writeEntryIfDifferent<word>(os, "p", "p", pName_);
    T0_.writeEntry("T0", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

makePatchTypeField
(
    fvPatchScalarField,
    isentropicTotalTemperatureFvPatchScalarField
);

} // End namespace Foam


// ************************************************************************* //
