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

#include "isentropicTotalPressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "basicThermo.H"
#include "isentropicTotalTemperatureFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::isentropicTotalPressureFvPatchScalarField::
isentropicTotalPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    UName_("U"),
    TName_("T"),
    p0_(p.size(), 0.0)
{}


Foam::isentropicTotalPressureFvPatchScalarField::
isentropicTotalPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    UName_(dict.lookupOrDefault<word>("U", "U")),
    TName_(dict.lookupOrDefault<word>("T", "T")),
    p0_("p0", dict, p.size())
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
        fvPatchField<scalar>::operator=(p0_);
    }
}


Foam::isentropicTotalPressureFvPatchScalarField::
isentropicTotalPressureFvPatchScalarField
(
    const isentropicTotalPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    UName_(ptf.UName_),
    TName_(ptf.TName_),
    p0_(ptf.p0_, mapper)
{}


Foam::isentropicTotalPressureFvPatchScalarField::
isentropicTotalPressureFvPatchScalarField
(
    const isentropicTotalPressureFvPatchScalarField& ptf
)
:
    fixedValueFvPatchScalarField(ptf),
    UName_(ptf.UName_),
    TName_(ptf.TName_),
    p0_(ptf.p0_)
{}


Foam::isentropicTotalPressureFvPatchScalarField::
isentropicTotalPressureFvPatchScalarField
(
    const isentropicTotalPressureFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    UName_(ptf.UName_),
    TName_(ptf.TName_),
    p0_(ptf.p0_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::isentropicTotalPressureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchScalarField::autoMap(m);
    p0_.autoMap(m);
}


void Foam::isentropicTotalPressureFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchScalarField::rmap(ptf, addr);

    const isentropicTotalPressureFvPatchScalarField& tiptf =
        refCast<const isentropicTotalPressureFvPatchScalarField>(ptf);

    p0_.rmap(tiptf.p0_, addr);
}


void Foam::isentropicTotalPressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Get velocity
    const fvPatchVectorField& U =
        patch().lookupPatchField<volVectorField, scalar>(UName_);

    // Get temperature
    const fvPatchScalarField& T =
        patch().lookupPatchField<volScalarField, scalar>(TName_);

    const basicThermo& thermo =
        db().lookupObject<basicThermo>("thermophysicalProperties");

    scalarField Cp = thermo.Cp(T, patch().index());
    scalarField Cv = thermo.Cv(T, patch().index());

    scalarField gamma = Cp/Cv;
    scalarField R = Cp - Cv;

    scalarField Ma = -(patch().nf() & U)/sqrt(gamma*R*T);

    scalarField a = 1 + 0.5*(gamma - 1)*sqr(Ma);

    operator==(p0_*pow(a, -gamma/(gamma - 1)));

    fixedValueFvPatchScalarField::updateCoeffs();
}

void Foam::isentropicTotalPressureFvPatchScalarField::updateCoeffs
(
    const vectorField& Up
)
{
    updateCoeffs();
}


void Foam::isentropicTotalPressureFvPatchScalarField::write
(
    Ostream& os
) const
{
    fixedValueFvPatchScalarField::write(os);
    writeEntryIfDifferent<word>(os, "T", "T", TName_);
    p0_.writeEntry("p0", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

makePatchTypeField
(
    fvPatchScalarField,
    isentropicTotalPressureFvPatchScalarField
);

} // End namespace Foam


// ************************************************************************* //
