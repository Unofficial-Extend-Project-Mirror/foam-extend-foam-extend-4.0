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

#include "greyDiffusiveRadiationMixedFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

#include "fvDOM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{
// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void
greyDiffusiveRadiationMixedFvPatchScalarField::calcSumOutgoingAngles() const
{
    if (sumOutgoingAnglesPtr_)
    {
        FatalErrorIn
        (
            "wideBandDiffusiveRadiationMixedFvPatchScalarField"
            "::calcSumOutgoingAngles()"
        )   << "sumOutgoingAnglesPtr_ already calculated"
            << abort(FatalError);
    }

    sumOutgoingAnglesPtr_ = new scalarField(this->size(), 0.);
    scalarField& sumOutgoingAngles = *sumOutgoingAnglesPtr_;


    // Get access to radiation model, and recast as fvDOM
    const radiationModel& radiation =
        db().lookupObject<radiationModel>("radiationProperties");

    const fvDOM& dom = dynamic_cast<const fvDOM&>(radiation);


    for(label rayI = 0; rayI < dom.nRay(); rayI++)
    {
        // Calculate cosine of angle between face and ray
        scalarField outgoingAngles = this->patch().nf() & dom.IRay(rayI).dAve();

        // For outgoing rays, outgoingAngles will be negative
        sumOutgoingAngles += neg(outgoingAngles)*(-outgoingAngles);
    }
}

const Foam::scalarField&
greyDiffusiveRadiationMixedFvPatchScalarField::sumOutgoingAngles() const
{
    if (!sumOutgoingAnglesPtr_)
    {
        calcSumOutgoingAngles();
    }

    return *sumOutgoingAnglesPtr_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

greyDiffusiveRadiationMixedFvPatchScalarField::
greyDiffusiveRadiationMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    TName_("undefinedT"),
    emissivity_(0.0),
    sumOutgoingAnglesPtr_(NULL)
{
    refValue() = 0.0;
    refGrad() = 0.0;
    valueFraction() = 1.0;
}


greyDiffusiveRadiationMixedFvPatchScalarField::
greyDiffusiveRadiationMixedFvPatchScalarField
(
    const greyDiffusiveRadiationMixedFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    TName_(ptf.TName_),
    emissivity_(ptf.emissivity_),
    sumOutgoingAnglesPtr_(NULL)
{}


greyDiffusiveRadiationMixedFvPatchScalarField::
greyDiffusiveRadiationMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    TName_(dict.lookup("T")),
    emissivity_(readScalar(dict.lookup("emissivity"))),
    sumOutgoingAnglesPtr_(NULL)
{
    if (dict.found("refValue"))
    {
        refValue() = scalarField("value", dict, p.size());

        fvPatchScalarField::operator=
        (
            refValue()
        );

        refGrad() = scalarField("refGradient", dict, p.size());
        valueFraction() = scalarField("valueFraction", dict, p.size());
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
        valueFraction() = 1;

        fvPatchScalarField::operator=(refValue());
    }
}


greyDiffusiveRadiationMixedFvPatchScalarField::
greyDiffusiveRadiationMixedFvPatchScalarField
(
    const greyDiffusiveRadiationMixedFvPatchScalarField& ptf
)
:
    mixedFvPatchScalarField(ptf),
    TName_(ptf.TName_),
    emissivity_(ptf.emissivity_),
    sumOutgoingAnglesPtr_(NULL)
{}


greyDiffusiveRadiationMixedFvPatchScalarField::
greyDiffusiveRadiationMixedFvPatchScalarField
(
    const greyDiffusiveRadiationMixedFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    TName_(ptf.TName_),
    emissivity_(ptf.emissivity_),
    sumOutgoingAnglesPtr_(NULL)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void greyDiffusiveRadiationMixedFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }




    const label patchI = this->patch().index();

    // Access radiation model
    const radiationModel& radiation =
        db().lookupObject<radiationModel>("radiationProperties");

    const fvDOM& dom = dynamic_cast<const fvDOM&>(radiation);

    if (dom.nLambda() == 0)
    {
        FatalErrorIn
        (
            ""
            "wideBandDiffusiveRadiationMixedFvPatchScalarField::updateCoeffs"
        )   << " a non-grey boundary condition is used with a grey "
            << "absorption model" << nl << exit(FatalError);
    }

    // Get rayId and lambda Id for this ray
    label rayId = -1;
    label lambdaId = -1;
    dom.setRayIdLambdaId(dimensionedInternalField().name(), rayId, lambdaId);

    // Make shortcut to ray belonging to this field
    const radiativeIntensityRay& ray = dom.IRay(rayId);

    // Access incoming radiation for this patch for this band
    const scalarField& Qin = dom.Qin(lambdaId)[patchI];

    // Access black body radiation for this patch for this band
    const scalarField& Eb =
        dom.blackBody().bLambda(lambdaId).boundaryField()[patchI];

    // Get face normals
    vectorField nHat = patch().nf();

    // Calculate cos of incoming angle of current ray with every face
    scalarField incomingAngle = this->patch().nf() & ray.dAve();

    // Set to zeroGradient (=0; incomingAngle > 0) for faces with incoming rays
    // and to fixedValue (=1; incomingAngle < 0) for outgoing rays
    this->valueFraction() = neg(incomingAngle);

    // Set intensity value for fixedValue part (see reference in header file)
    this->refValue() =
        emissivity_*Eb + ((1 - emissivity_)*Qin)/sumOutgoingAngles();

    // Update boundary field now, so values for incoming and outgoing rays
    // are in balance
    scalarField::operator=
    (
        this->valueFraction()*this->refValue()
      +
        (1.0 - this->valueFraction())*
        (
            this->patchInternalField()
          + this->refGrad()/this->patch().deltaCoeffs()
        )
    );

    mixedFvPatchScalarField::updateCoeffs();
}


void greyDiffusiveRadiationMixedFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    os.writeKeyword("T") << TName_ << token::END_STATEMENT << nl;
    os.writeKeyword("emissivity") << emissivity_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    greyDiffusiveRadiationMixedFvPatchScalarField
);

} // End namespace radiation
} // End namespace Foam

// ************************************************************************* //
