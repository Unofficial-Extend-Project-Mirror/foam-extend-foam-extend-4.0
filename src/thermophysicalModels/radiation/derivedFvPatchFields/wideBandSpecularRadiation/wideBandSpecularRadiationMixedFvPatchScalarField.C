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

#include "radiationModel.H"
#include "fvDOM.H"
#include "wideBandSpecularRadiationMixedFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{
// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void wideBandSpecularRadiationMixedFvPatchScalarField::
calcReceivedRayIDs() const
{
    if (receivedRayIDPtr_)
    {
        FatalErrorIn
        (
            "wideBandSpecularRadiationMixedFvPatchScalarField"
            "::calcreceivedRayIDs()"
        )   << "receivedRayIDPtr already calculated"
            << abort(FatalError);
    }

    receivedRayIDPtr_ = new labelListList(this->size());
    labelListList& receivedRayID = *receivedRayIDPtr_;

    // Get access to radiation model, and recast as fvDOM
    const radiationModel& radiation =
        this->db().lookupObject<radiationModel>("radiationProperties");

    const fvDOM& dom(refCast<const fvDOM>(radiation));

    // Get rayId and lambda Id for this ray
    label rayId = -1;
    label lambdaId = -1;
    dom.setRayIdLambdaId(dimensionedInternalField().name(), rayId, lambdaId);

    // Get index of active ray
    word rayAndBand = this->dimensionedInternalField().name();
    rayAndBand =
        rayAndBand.substr(rayAndBand.find("_")+1, rayAndBand.size()-1);

    // Get all ray directions
    List<vector> dAve(dom.nRay());
    forAll (dAve, rayI)
    {
        dAve[rayI] = dom.IRay(rayI).dAve()/mag(dom.IRay(rayI).dAve());
    }

    // Get face normal vectors
    vectorField nHat = this->patch().nf();


    // For each face, and for each reflected ray, try to find another
    // ray that is better suited to pick it up.
    forAll(receivedRayID, faceI)
    {
        // Access the list of received rays for this face
        labelList& receivedRayIDs = receivedRayID[faceI];

        // Check whether ray is going into surface
        // -> no reflection
        if ((dAve[rayId] & nHat[faceI]) > 0)
        {
            receivedRayIDs.setSize(0);
            continue;
        }


        // For each face, initialize list of picked up
        // rays to maximum possible size
        receivedRayIDs.setSize(dom.nRay() - 1);

        // Count actual number of receiving rays
        label nReceiving = 0;

        forAll(dAve, incomingRayI)
        {
            // Calculate reflected ray direction for rayI
            vector dReflected =
                dAve[incomingRayI]
              - 2*(dAve[incomingRayI] & nHat[faceI])*nHat[faceI];

            // Get dot product with this ray
            scalar dotProductThisRay = dAve[rayId] & dReflected;

            // Assume this ray is closest to the reflected ray
            bool closest = true;

            // And look through all other rays to find a better suited one
            forAll(dAve, receivingRayI)
            {
                if (receivingRayI == rayId)
                {
                    // Do not compare this ray with itself
                    continue;
                }

                scalar dotProductOtherRay = dAve[receivingRayI] & dReflected;
                if (dotProductThisRay < dotProductOtherRay)
                {
                    // If another ray is closer, stop search
                    closest = false;
                    break;
                }
            }

            if (closest)
            {
                // Could not find better suited ray, so this ray needs to
                // pick it up. Add incoming ray to list for this face
                receivedRayIDs[nReceiving] = incomingRayI;
                nReceiving ++;
            }

        }

        // Resize list of picked up rays for this face
        receivedRayIDs.setSize(nReceiving);
    }

    // Sanity check on last ray
    if (rayId == (dom.nRay() - 1))
    {

        forAll(nHat, faceI)
        {
            label incomingRays = 0;
            label pickedUpRays = 0;

            for (label rayI = 0; rayI < dom.nRay(); rayI++)
            {
                // Is this ray going into the wall?
                if ((dAve[rayI] &  nHat[faceI]) > 0)
                {
                    incomingRays++;
                }

                // How many rays are picked up by this ray?
                const wideBandSpecularRadiationMixedFvPatchScalarField& rayBC =
                    refCast
                    <
                        const wideBandSpecularRadiationMixedFvPatchScalarField
                    >
                    (
                        dom.IRayLambda
                        (
                            rayI,
                            lambdaId
                        ).boundaryField()[patch().index()]
                    );

                const labelListList& receivedRaysList = rayBC.receivedRayIDs();

                pickedUpRays += receivedRaysList[faceI].size();
            }

            if (incomingRays != pickedUpRays)
            {
                FatalErrorIn
                (
                    "wideBandSpecularRadiationMixedFvPatchScalarField"
                    "::calcreceivedRayIDs()"
                )   << "Sanity checked failed on patch " << patch().name()
                    << " face " << faceI << nl
                    << "number of rays with direction into wall: "
                    << incomingRays
                    << ", number of reflected rays picked up by outgoing rays: "
                    << pickedUpRays << nl
                    << abort(FatalError);
            }
        }
    }

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

wideBandSpecularRadiationMixedFvPatchScalarField::
wideBandSpecularRadiationMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchField<scalar>(p, iF),
    receivedRayIDPtr_(NULL)
{
    this->refValue() = 0;
    this->refGrad() = 0;
    this->valueFraction() = 0;
}


wideBandSpecularRadiationMixedFvPatchScalarField::
wideBandSpecularRadiationMixedFvPatchScalarField
(
    const wideBandSpecularRadiationMixedFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<scalar>(ptf, p, iF, mapper),
    receivedRayIDPtr_(NULL)
{}


wideBandSpecularRadiationMixedFvPatchScalarField::
wideBandSpecularRadiationMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<scalar>(p, iF),
    receivedRayIDPtr_(NULL)
{
    this->refValue() = 0;
    this->refGrad() = 0;
    this->valueFraction() = 0;

    if (dict.found("patchType"))
    {
        word& pType = const_cast<word &>(this->patchType());
        pType = word(dict.lookup("patchType"));
}

    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        fvPatchScalarField::operator=(this->refValue());
    }
}


wideBandSpecularRadiationMixedFvPatchScalarField::
wideBandSpecularRadiationMixedFvPatchScalarField
(
    const wideBandSpecularRadiationMixedFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchField<scalar>(ptf, iF),
    receivedRayIDPtr_(NULL)
{}


wideBandSpecularRadiationMixedFvPatchScalarField::
wideBandSpecularRadiationMixedFvPatchScalarField
(
    const wideBandSpecularRadiationMixedFvPatchScalarField& ptf
)
:
    mixedFvPatchField<scalar>(ptf),
    receivedRayIDPtr_(NULL)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::labelListList&
wideBandSpecularRadiationMixedFvPatchScalarField::receivedRayIDs() const
{
    if (!receivedRayIDPtr_)
    {
        calcReceivedRayIDs();
    }

    return *receivedRayIDPtr_;
}


// Evaluate the field on the patch
void wideBandSpecularRadiationMixedFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    vectorField nHat = this->patch().nf();

    // Get access to radiation model, and recast as fvDOM
    const radiationModel& radiation =
        this->db().lookupObject<radiationModel>("radiationProperties");

    const fvDOM& dom(refCast<const fvDOM>(radiation));

    // Get rayId and lambda Id for this ray
    label rayId = -1;
    label lambdaId = -1;
    dom.setRayIdLambdaId(dimensionedInternalField().name(), rayId, lambdaId);

    const labelListList& receivedRayIDs = this->receivedRayIDs();

    scalarField IValue(this->size(), 0);

    labelList faceCellIDs = this->patch().faceCells();

    // Loop over all faces and add values from all received rays
    forAll(receivedRayIDs, faceI)
    {
        if (receivedRayIDs[faceI].size() == 0)
        {
            // Ray goes into face -> act as zeroGradient
            this->valueFraction()[faceI] = 0;
        }
        else
        {
            // Ray goes out of face -> act as fixedValue
            this->valueFraction()[faceI] = 1;

            // Pick up all reflected rays
            forAll(receivedRayIDs[faceI], receivedRayI)
            {
                // Get ray from object registry
                IValue[faceI] +=
                    dom.IRayLambda
                    (
                        receivedRayIDs[faceI][receivedRayI],
                        lambdaId
                    ).internalField()[faceCellIDs[faceI]];
            }
        }
    }

    // Set value for patch
    this->refValue() = IValue;

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

    mixedFvPatchField<scalar>::updateCoeffs();
}


void wideBandSpecularRadiationMixedFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void wideBandSpecularRadiationMixedFvPatchScalarField::operator=
(
    const fvPatchScalarField& ptf
)
{
    fvPatchScalarField::operator=
    (
        this->valueFraction()*this->refValue()
        + (1 - this->valueFraction())*ptf
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
   fvPatchScalarField,
   wideBandSpecularRadiationMixedFvPatchScalarField
);


} // End namespace radiation
} // End namespace Foam

// ************************************************************************* //
