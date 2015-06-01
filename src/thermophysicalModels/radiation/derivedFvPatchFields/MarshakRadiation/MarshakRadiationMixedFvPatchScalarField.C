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

#include "MarshakRadiationMixedFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

#include "fvc.H"
#include "radiationModel.H"
#include "radiationConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::MarshakRadiationFvPatchScalarField::MarshakRadiationFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    TName_("undefined"),
    emissivity_(0.0)
{
    refValue() = 0.0;
    refGrad() = 0.0;
    valueFraction() = 0.0;
}


Foam::MarshakRadiationFvPatchScalarField::MarshakRadiationFvPatchScalarField
(
    const MarshakRadiationFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    TName_(ptf.TName_),
    emissivity_(ptf.emissivity_)
{}


Foam::MarshakRadiationFvPatchScalarField::MarshakRadiationFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    TName_(dict.lookup("T")),
    emissivity_(readScalar(dict.lookup("emissivity")))
{
    if (dict.found("value"))
    {
        refValue() = scalarField("value", dict, p.size());

        fvPatchScalarField::operator=
        (
            refValue()
        );

        refGrad() = 0;
        valueFraction() = 1;
    }
    else
    {
        // No value given. Restart as fixedValue b.c.

        // Bugfix: Do not initialize from temperautre because it is unavailable
        // when running, e.g. decomposePar and loading radiation as
        // shared library. Initialize to zero instead.
        // 26 Mar 2014 - DC
        refValue() = 0;

        refGrad() = 0;
        valueFraction() = 0;

        fvPatchScalarField::operator=(refValue());
    }
}


Foam::MarshakRadiationFvPatchScalarField::MarshakRadiationFvPatchScalarField
(
    const MarshakRadiationFvPatchScalarField& ptf
)
:
    mixedFvPatchScalarField(ptf),
    TName_(ptf.TName_),
    emissivity_(ptf.emissivity_)
{}


Foam::MarshakRadiationFvPatchScalarField::MarshakRadiationFvPatchScalarField
(
    const MarshakRadiationFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    TName_(ptf.TName_),
    emissivity_(ptf.emissivity_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::MarshakRadiationFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    scalarField::autoMap(m);
}


void Foam::MarshakRadiationFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    mixedFvPatchScalarField::rmap(ptf, addr);
}


void Foam::MarshakRadiationFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Temperature field
    const scalarField& Tp =
        lookupPatchField<volScalarField, scalar>(TName_);

    // Re-calc reference value
    refValue() = 4.0*radiation::sigmaSB.value()*pow4(Tp);

    // Diffusion coefficient - created by radiation model's ::updateCoeffs()
    const scalarField& gamma =
        lookupPatchField<volScalarField, scalar>("gammaRad");

    const scalar Ep = emissivity_/(2.0*(2.0 - emissivity_));

    // Set value fraction
    valueFraction() = 1.0/(1.0 + gamma*patch().deltaCoeffs()/Ep);

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::MarshakRadiationFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("T") << TName_ << token::END_STATEMENT << nl;
    os.writeKeyword("emissivity") << emissivity_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        MarshakRadiationFvPatchScalarField
    );
}


// ************************************************************************* //
