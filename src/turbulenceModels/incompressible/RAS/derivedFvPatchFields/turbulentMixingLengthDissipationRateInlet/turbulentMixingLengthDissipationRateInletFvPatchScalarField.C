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

#include "turbulentMixingLengthDissipationRateInletFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "RASModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

turbulentMixingLengthDissipationRateInletFvPatchScalarField::
turbulentMixingLengthDissipationRateInletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    inletOutletFvPatchScalarField(p, iF),
    mixingLength_(0.0),
    phiName_("undefined-phi"),
    kName_("undefined-k")
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 0.0;
}

turbulentMixingLengthDissipationRateInletFvPatchScalarField::
turbulentMixingLengthDissipationRateInletFvPatchScalarField
(
    const turbulentMixingLengthDissipationRateInletFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    inletOutletFvPatchScalarField(ptf, p, iF, mapper),
    mixingLength_(ptf.mixingLength_),
    phiName_(ptf.phiName_),
    kName_(ptf.kName_)
{}

turbulentMixingLengthDissipationRateInletFvPatchScalarField::
turbulentMixingLengthDissipationRateInletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    inletOutletFvPatchScalarField(p, iF),
    mixingLength_(readScalar(dict.lookup("mixingLength"))),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    kName_(dict.lookupOrDefault<word>("k", "k"))
{
    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 0.0;
}

turbulentMixingLengthDissipationRateInletFvPatchScalarField::
turbulentMixingLengthDissipationRateInletFvPatchScalarField
(
    const turbulentMixingLengthDissipationRateInletFvPatchScalarField& ptf
)
:
    inletOutletFvPatchScalarField(ptf),
    mixingLength_(ptf.mixingLength_),
    phiName_(ptf.phiName_),
    kName_(ptf.kName_)
{}

turbulentMixingLengthDissipationRateInletFvPatchScalarField::
turbulentMixingLengthDissipationRateInletFvPatchScalarField
(
    const turbulentMixingLengthDissipationRateInletFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    inletOutletFvPatchScalarField(ptf, iF),
    mixingLength_(ptf.mixingLength_),
    phiName_(ptf.phiName_),
    kName_(ptf.kName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void turbulentMixingLengthDissipationRateInletFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Lookup Cmu corresponding to the turbulence model selected
    const RASModel& rasModel = db().lookupObject<RASModel>("RASProperties");

    const scalar Cmu =
        rasModel.coeffDict().lookupOrDefault<scalar>("Cmu", 0.09);

    const scalar Cmu75 = pow(Cmu, 0.75);

    const fvPatchScalarField& kp =
        lookupPatchField<volScalarField, scalar>(kName_);

    const fvsPatchScalarField& phip =
        lookupPatchField<surfaceScalarField, scalar>(phiName_);

    this->refValue() = Cmu75*kp*sqrt(kp)/mixingLength_;
    this->valueFraction() = neg(phip);

    inletOutletFvPatchScalarField::updateCoeffs();
}


void turbulentMixingLengthDissipationRateInletFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    os.writeKeyword("mixingLength")
        << mixingLength_ << token::END_STATEMENT << nl;
    os.writeKeyword("phi") << phiName_ << token::END_STATEMENT << nl;
    os.writeKeyword("k") << kName_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    turbulentMixingLengthDissipationRateInletFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
