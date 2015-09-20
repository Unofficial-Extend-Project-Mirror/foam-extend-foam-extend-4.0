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
    fixedGradientFvPatchScalarField(p, iF),
    q_(p.size(), 0.0),
    KappaName_("undefined-Kappa")
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
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    q_(ptf.q_, mapper),
    KappaName_(ptf.KappaName_)
{}


Foam::solidWallHeatFluxTemperatureFvPatchScalarField::
solidWallHeatFluxTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF, dict),
    q_("q", dict, p.size()),
    KappaName_(dict.lookup("Kappa"))
{}


Foam::solidWallHeatFluxTemperatureFvPatchScalarField::
solidWallHeatFluxTemperatureFvPatchScalarField
(
    const solidWallHeatFluxTemperatureFvPatchScalarField& tppsf
)
:
    fixedGradientFvPatchScalarField(tppsf),
    q_(tppsf.q_),
    KappaName_(tppsf.KappaName_)
{}


Foam::solidWallHeatFluxTemperatureFvPatchScalarField::
solidWallHeatFluxTemperatureFvPatchScalarField
(
    const solidWallHeatFluxTemperatureFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(tppsf, iF),
    q_(tppsf.q_),
    KappaName_(tppsf.KappaName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::solidWallHeatFluxTemperatureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedGradientFvPatchScalarField::autoMap(m);
    q_.autoMap(m);
}


void Foam::solidWallHeatFluxTemperatureFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedGradientFvPatchScalarField::rmap(ptf, addr);

    const solidWallHeatFluxTemperatureFvPatchScalarField& hfptf =
        refCast<const solidWallHeatFluxTemperatureFvPatchScalarField>(ptf);

    q_.rmap(hfptf.q_, addr);
}


Foam::tmp<Foam::scalarField>
Foam::solidWallHeatFluxTemperatureFvPatchScalarField::Kappa() const
{
    const fvMesh& mesh = patch().boundaryMesh().mesh();

    if (mesh.objectRegistry::foundObject<volScalarField>(KappaName_))
    {
        return lookupPatchField<volScalarField, scalar>(KappaName_);
    }
    else if (mesh.objectRegistry::foundObject<volSymmTensorField>(KappaName_))
    {
        const symmTensorField& KappaWall =
            lookupPatchField<volSymmTensorField, scalar>(KappaName_);

        vectorField n = patch().nf();

        return n & KappaWall & n;
    }
    else
    {
        FatalErrorIn
        (
            "solidWallHeatFluxTemperatureFvPatchScalarField::Kappa()"
            " const"
        )   << "Did not find field " << KappaName_
            << " on mesh " << mesh.name() << " patch " << patch().name()
            << endl
            << "Please set 'Kappa' to a valid volScalarField"
            << " or a valid volSymmTensorField." << exit(FatalError);

        return scalarField(0);
    }
}


void Foam::solidWallHeatFluxTemperatureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    gradient() = q_/Kappa();

    fixedGradientFvPatchScalarField::updateCoeffs();

    if (debug)
    {
        scalar Q = gSum(Kappa()*patch().magSf()*snGrad());

        Info<< patch().boundaryMesh().mesh().name() << ':'
            << patch().name() << ':'
            << this->dimensionedInternalField().name() << " :"
            << " heatFlux:" << Q
            << " walltemperature "
            << " min:" << gMin(*this)
            << " max:" << gMax(*this)
            << " avg:" << gAverage(*this)
            << endl;
    }
}


void Foam::solidWallHeatFluxTemperatureFvPatchScalarField::write
(
    Ostream& os
) const
{
    fixedGradientFvPatchScalarField::write(os);
    q_.writeEntry("q", os);
    os.writeKeyword("Kappa") << KappaName_ << token::END_STATEMENT << nl;
    this->writeEntry("value", os);
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
