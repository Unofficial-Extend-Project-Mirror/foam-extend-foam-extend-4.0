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

#include "solidWallHeatFluxTemperatureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidWallHeatFluxTemperatureFvPatchScalarField::
solidWallHeatFluxTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    q_(p.size(), 0.0),
    KName_("undefined-K")
{}


Foam::solidWallHeatFluxTemperatureFvPatchScalarField::
solidWallHeatFluxTemperatureFvPatchScalarField
(
    const solidWallHeatFluxTemperatureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    q_(ptf.q_, mapper),
    KName_(ptf.KName_)
{}


Foam::solidWallHeatFluxTemperatureFvPatchScalarField::
solidWallHeatFluxTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    q_("q", dict, p.size()),
    KName_(dict.lookup("K"))
{}


Foam::solidWallHeatFluxTemperatureFvPatchScalarField::
solidWallHeatFluxTemperatureFvPatchScalarField
(
    const solidWallHeatFluxTemperatureFvPatchScalarField& tppsf
)
:
    fixedValueFvPatchScalarField(tppsf),
    q_(tppsf.q_),
    KName_(tppsf.KName_)
{}


Foam::solidWallHeatFluxTemperatureFvPatchScalarField::
solidWallHeatFluxTemperatureFvPatchScalarField
(
    const solidWallHeatFluxTemperatureFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(tppsf, iF),
    q_(tppsf.q_),
    KName_(tppsf.KName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::solidWallHeatFluxTemperatureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchScalarField::autoMap(m);
    q_.autoMap(m);
}


void Foam::solidWallHeatFluxTemperatureFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchScalarField::rmap(ptf, addr);

    const solidWallHeatFluxTemperatureFvPatchScalarField& hfptf =
        refCast<const solidWallHeatFluxTemperatureFvPatchScalarField>(ptf);

    q_.rmap(hfptf.q_, addr);
}


void Foam::solidWallHeatFluxTemperatureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const scalarField& Kw =
        patch().lookupPatchField<volScalarField, scalar>(KName_);

    const fvPatchScalarField& Tw = *this;

    operator==(q_/(patch().deltaCoeffs()*Kw) + Tw.patchInternalField());

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::solidWallHeatFluxTemperatureFvPatchScalarField::write
(
    Ostream& os
) const
{
    fixedValueFvPatchScalarField::write(os);
    q_.writeEntry("q", os);
    os.writeKeyword("K") << KName_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        solidWallHeatFluxTemperatureFvPatchScalarField
    );
}

// ************************************************************************* //
