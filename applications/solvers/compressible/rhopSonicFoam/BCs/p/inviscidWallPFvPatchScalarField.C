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

#include "inviscidWallPFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inviscidWallPFvPatchScalarField::inviscidWallPFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    fluxFraction_(1)
{}


inviscidWallPFvPatchScalarField::inviscidWallPFvPatchScalarField
(
    const inviscidWallPFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    fluxFraction_(ptf.fluxFraction_)
{}


inviscidWallPFvPatchScalarField::inviscidWallPFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),
    fluxFraction_(readScalar(dict.lookup("fluxFraction")))
{
    if (dict.found("gradient"))
    {
        gradient() = scalarField("gradient", dict, p.size());
        fixedGradientFvPatchScalarField::updateCoeffs();
        fixedGradientFvPatchScalarField::evaluate();
    }
    else
    {
        fvPatchField<scalar>::operator=(patchInternalField());
        gradient() = 0.0;
    }

    if (fluxFraction_<0.0 || fluxFraction_ > 1.0)
    {
        FatalIOErrorIn
        (
            "inviscidWallPFvPatchScalarField::"
            "supersonicFreeStreamFvPatchVectorField"
            "(const fvPatch&, const scalarField&, const dictionary&)",
            dict
        )   << "    unphysical fluxFraction specified (< 0.0 or > 1.0)"
            << exit(FatalIOError);
    }

}


inviscidWallPFvPatchScalarField::inviscidWallPFvPatchScalarField
(
    const inviscidWallPFvPatchScalarField& wbppsf
)
:
    fixedGradientFvPatchScalarField(wbppsf),
    fluxFraction_(wbppsf.fluxFraction_)
{}


inviscidWallPFvPatchScalarField::inviscidWallPFvPatchScalarField
(
    const inviscidWallPFvPatchScalarField& wbppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(wbppsf, iF),
    fluxFraction_(wbppsf.fluxFraction_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void inviscidWallPFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvPatchField<vector>& rhoUp =
        lookupPatchField<volVectorField, vector>("rhoU");

    const fvsPatchField<scalar>& phip =
        lookupPatchField<surfaceScalarField, scalar>("phi");

    const fvsPatchField<scalar>& rAp =
        lookupPatchField<surfaceScalarField, scalar>("rrhoUAf");

    gradient() = (fluxFraction_*phip - (patch().Sf() & rhoUp))/
                 (rAp*patch().magSf());

    fixedGradientFvPatchScalarField::updateCoeffs();
}


void inviscidWallPFvPatchScalarField::write(Ostream& os) const
{
    fixedGradientFvPatchScalarField::write(os);
    os.writeKeyword("fluxFraction")
        << fluxFraction_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, inviscidWallPFvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
