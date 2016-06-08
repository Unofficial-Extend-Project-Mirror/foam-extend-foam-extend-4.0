/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2007 Hrvoje Jasak
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Class
    cohesiveZoneFvPatchVectorField

Description

\*---------------------------------------------------------------------------*/

#include "cohesiveZoneFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "transformField.H"
#include "volFields.H"
#include "pointFields.H"
#include "constitutiveModel.H"
// #include "solidSolver.H"
#include "wallFvPatch.H"
#include "componentMixedPointPatchFields.H"

#include "tractionBoundaryGradient.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

cohesiveZoneFvPatchVectorField
::cohesiveZoneFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    directionMixedDisplacementFvPatchVectorField(p, iF),
//     totRefValue_(p.size(), vector::zero),
//     totalFieldName_("undefined"),
    cohesiveLawPtr_(NULL),
    crackIndicator_(p.size(), 0.0),
    crazeIndicator_(p.size(), 0.0),
    relaxationFactor_(1.0),
    separationDistance_(p.size(), vector::zero),
    oldSeparationDistance_(p.size(), vector::zero),
    unloadingSeparationDistance_(p.size(), vector::zero),
    explicitSeparationDistance_(false),
    curTimeIndex_(-1),
    traction_(p.size(), vector::zero),
    initiationTraction_(p.size(), vector::zero),
    beta_(0.0),
    breakOnlyOneFace_(false),
    timeSeries_()
{}


cohesiveZoneFvPatchVectorField
::cohesiveZoneFvPatchVectorField
(
    const cohesiveZoneFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    directionMixedDisplacementFvPatchVectorField(ptf, p, iF, mapper),
    totRefValue_(ptf.totRefValue_),
//     totalFieldName_(ptf.totalFieldName_),
    cohesiveLawPtr_(ptf.cohesiveLawPtr_),
    crackIndicator_(ptf.crackIndicator_),
    crazeIndicator_(ptf.crazeIndicator_),
    relaxationFactor_(ptf.relaxationFactor_),
    separationDistance_(ptf.separationDistance_),
    oldSeparationDistance_(ptf.oldSeparationDistance_),
    unloadingSeparationDistance_(ptf.unloadingSeparationDistance_),
    explicitSeparationDistance_(ptf.explicitSeparationDistance_),
    curTimeIndex_(ptf.curTimeIndex_),
    traction_(ptf.initiationTraction_),
    initiationTraction_(ptf.initiationTraction_),
    beta_(ptf.beta_),
    breakOnlyOneFace_(ptf.breakOnlyOneFace_),
    timeSeries_()
{
    if (ptf.timeSeries_.valid())
    {
        timeSeries_.set(new interpolationTable<vector>(ptf.timeSeries_()));
    }
}


cohesiveZoneFvPatchVectorField
::cohesiveZoneFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    directionMixedDisplacementFvPatchVectorField(p, iF),
    totRefValue_(p.size(), vector::zero),
//     totalFieldName_(dict.lookup("totalFieldName")),
    cohesiveLawPtr_
    (
        simpleCohesiveLaw::New(dict.lookup("simpleCohesiveLaw"), dict).ptr()
    ),
    crackIndicator_(p.size(), 0.0),
    crazeIndicator_(p.size(), 0.0),
    relaxationFactor_(readScalar(dict.lookup("relaxationFactor"))),
    separationDistance_(p.size(), vector::zero),
    oldSeparationDistance_(p.size(), vector::zero),
    unloadingSeparationDistance_(p.size(), vector::zero),
    explicitSeparationDistance_(dict.lookup("explicitSeparationDistance")),
    curTimeIndex_(-1),
    traction_(p.size(), vector::zero),
    initiationTraction_(p.size(), vector::zero),
    beta_(readScalar(dict.lookup("beta"))),
    breakOnlyOneFace_(dict.lookup("breakOnlyOneFace")),
    timeSeries_()
{
    Info << "Creating " << cohesiveZoneFvPatchVectorField::typeName
        << " boundary" << endl;

    if (dict.found("totRefValue"))
    {
        totRefValue_ = vectorField("totRefValue", dict, p.size());
    }
    else
    {
        totRefValue_ = vector::zero;
    }

    if (dict.found("refValue"))
    {
        this->refValue() = vectorField("refValue", dict, p.size());
    }
    else
    {
        this->refValue() = vector::zero;
    }

    if (dict.found("refGradient"))
    {
        this->refGrad() = vectorField("refGradient", dict, p.size());
    }
    else
    {
        this->refGrad() = vector::zero;
    }

    if (dict.found("valueFraction"))
    {
        this->valueFraction() =
            symmTensorField("valueFraction", dict, p.size());
    }
    else
    {
        if (patch().type() == wallFvPatch::typeName)
        {
            this->valueFraction() = I;
        }
        else
        {
            // Symmetry plane
            vectorField n = this->patch().nf();
            this->valueFraction() = sqr(n);
        }
    }

    if (dict.found("value"))
    {
        Field<vector>::operator=(vectorField("value", dict, p.size()));
    }
    else
    {
        directionMixedFvPatchVectorField::evaluate();

//         Field<vector> normalValue = transform(valueFraction(), refValue());

//         Field<vector> gradValue =
//             this->patchInternalField() + refGrad()/this->patch().deltaCoeffs();

//         Field<vector> transformGradValue =
//             transform(I - valueFraction(), gradValue);

//         Field<vector>::operator=(normalValue + transformGradValue);
    }

    if (dict.found("crackIndicator"))
    {
        crackIndicator_ = scalarField("crackIndicator", dict, p.size());
    }

    if (dict.found("crazeIndicator"))
    {
        crazeIndicator_ = scalarField("crazeIndicator", dict, p.size());
    }

    if (dict.found("separationDistance"))
    {
        separationDistance_ =
            vectorField("separationDistance", dict, p.size());
    }

    if (dict.found("oldSeparationDistance"))
    {
        oldSeparationDistance_ =
            vectorField("oldSeparationDistance", dict, p.size());
    }

    if (dict.found("unloadingSeparationDistance"))
    {
        unloadingSeparationDistance_ =
            vectorField("unloadingSeparationDistance", dict, p.size());
    }

    if (dict.found("traction"))
    {
        traction_ =
            vectorField("traction", dict, p.size());
    }

    if (dict.found("initiationTraction"))
    {
        initiationTraction_ =
            vectorField("initiationTraction", dict, p.size());
    }

    // Check if time varying
    if (dict.found("fileName") && dict.found("outOfBounds"))
    {
        timeSeries_.set(new interpolationTable<vector>(dict));
    }
}


cohesiveZoneFvPatchVectorField
::cohesiveZoneFvPatchVectorField
(
    const cohesiveZoneFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    directionMixedDisplacementFvPatchVectorField(ptf, iF),
    totRefValue_(ptf.totRefValue_),
//     totalFieldName_(ptf.totalFieldName_),
    cohesiveLawPtr_(ptf.cohesiveLawPtr_),
    crackIndicator_(ptf.crackIndicator_),
    crazeIndicator_(ptf.crazeIndicator_),
    relaxationFactor_(ptf.relaxationFactor_),
    separationDistance_(ptf.separationDistance_),
    oldSeparationDistance_(ptf.oldSeparationDistance_),
    unloadingSeparationDistance_(ptf.unloadingSeparationDistance_),
    explicitSeparationDistance_(ptf.explicitSeparationDistance_),
    curTimeIndex_(ptf.curTimeIndex_),
    traction_(ptf.traction_),
    initiationTraction_(ptf.initiationTraction_),
    beta_(ptf.beta_),
    breakOnlyOneFace_(ptf.breakOnlyOneFace_),
    timeSeries_()
{
    if (ptf.timeSeries_.valid())
    {
        timeSeries_.set(new interpolationTable<vector>(ptf.timeSeries_()));
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Map from self
void cohesiveZoneFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    if (cohesiveLawPtr_ == NULL)
    {
        FatalErrorIn("cohesiveZoneFvPatchVectorField::autoMap")
            << "NULL cohesive law"
            << abort(FatalError);
    }

    directionMixedDisplacementFvPatchVectorField::autoMap(m);
    totRefValue_.autoMap(m);
    crackIndicator_.autoMap(m);
    crazeIndicator_.autoMap(m);
    separationDistance_.autoMap(m);
    oldSeparationDistance_.autoMap(m);
    unloadingSeparationDistance_.autoMap(m);
    traction_.autoMap(m);
    initiationTraction_.autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void cohesiveZoneFvPatchVectorField::rmap
(
    const fvPatchField<vector>& ptf,
    const labelList& addr
)
{
    directionMixedDisplacementFvPatchVectorField::rmap(ptf, addr);

    const cohesiveZoneFvPatchVectorField& dmptf =
        refCast<const cohesiveZoneFvPatchVectorField>(ptf);

    // No need to grab the cohesive zone pointer more than once
    if (!cohesiveLawPtr_)
    {
        cohesiveLawPtr_ = dmptf.cohesiveLawPtr_->clone().ptr();
    }

    totRefValue_.rmap(dmptf.totRefValue_, addr);
    crackIndicator_.rmap(dmptf.crackIndicator_, addr);
    crazeIndicator_.rmap(dmptf.crazeIndicator_, addr);
    relaxationFactor_ = dmptf.relaxationFactor_;
    separationDistance_.rmap(dmptf.separationDistance_, addr);
    oldSeparationDistance_.rmap(dmptf.oldSeparationDistance_, addr);
    unloadingSeparationDistance_.rmap
    (
        dmptf.unloadingSeparationDistance_, addr
    );
    explicitSeparationDistance_ = dmptf.explicitSeparationDistance_;
    curTimeIndex_ = dmptf.curTimeIndex_;
    traction_.rmap(dmptf.traction_, addr);
    initiationTraction_.rmap(dmptf.initiationTraction_, addr);
    breakOnlyOneFace_ = dmptf.breakOnlyOneFace_;

    if (dmptf.timeSeries_.valid())
    {
        timeSeries_.set(new interpolationTable<vector>(dmptf.timeSeries_()));
    }
}

tmp<scalarField> cohesiveZoneFvPatchVectorField::normalTraction() const
{
    vectorField n = patch().nf();
    tmp<scalarField> tNormalTraction
    (
        new scalarField(size(), 0.0)
    );

    scalarField nt = (n&traction_);

    for (label i=0; i<size(); i++)
    {
        tNormalTraction()[i] = nt[i];
    }

    return tNormalTraction;
}

tmp<scalarField> cohesiveZoneFvPatchVectorField::effectiveTraction
(
    const vectorField& traction,
    const vectorField& normal
) const
{
    scalarField normalTraction = (normal&traction);
    normalTraction *= pos(normalTraction);

    tmp<scalarField> effTraction
    (
        new scalarField(traction.size(), 0)
    );

    if (beta_ > SMALL)
    {
        effTraction() =
            sqrt
            (
                magSqr((I - normal*normal)&traction)/sqr(beta_)
              + sqr(normalTraction)
            );
    }
    else if (beta_ < -SMALL)
    {
        effTraction =
            mag((I - normal*normal) & traction);
    }
    else
    {
        effTraction() = normalTraction;
    }

    return effTraction;
}

scalar cohesiveZoneFvPatchVectorField::effectiveTraction
(
    const vector& traction,
    const vector& normal
) const
{
    scalar normalTraction = (normal & traction);
    normalTraction *= pos(normalTraction);

    scalar effTraction = 0;

    if (beta_ > SMALL)
    {
        effTraction =
            sqrt
            (
                magSqr((I - normal*normal)&traction)/sqr(beta_)
              + sqr(normalTraction)
            );
    }
    else if (beta_ < -SMALL)
    {
        effTraction =
            mag((I - normal*normal) & traction);
    }
    else
    {
        effTraction = normalTraction;
    }

    return effTraction;
}

label cohesiveZoneFvPatchVectorField::updateCrack()
{
    word DName = this->dimensionedInternalField().name();

    // Current stress
    const fvsPatchField<symmTensor>& sigma =
        patch().lookupPatchField<surfaceSymmTensorField, symmTensor>
        (
            "sigmaf"
        );

    // Patch stress
//     symmTensorField curSigma = sigma;
//     symmTensorField curSigma = sigma + DSigma;

    // Patch normal
    vectorField n = this->patch().nf();

    // Current traction
    traction_ = (n&sigma);

    // Current normal traction
    scalarField curNormalTraction = (n & (n & sigma));

    // Current tangential traction
    vectorField curTangentialTraction = ((I-n*n) & (n & sigma));

    scalarField effTraction = effectiveTraction(traction_, n);

    Info << "Max effective traction: "
        << gMax(effectiveTraction(traction_, n))
        << " (" << law().sigmaMax().value() << ")" << endl;

    // Select potential faces to breake
    SLList<label> facesToBreakList;
    SLList<scalar> facesToBreakTractionList;

    forAll(effTraction, faceI)
    {
        if
        (
            (magSqr(this->valueFraction()[faceI]) > 1-SMALL)
         && (effTraction[faceI] >= law().sigmaMax().value())
        )
        {
            facesToBreakList.insert(faceI);
            facesToBreakTractionList.insert(effTraction[faceI]);
        }
    }

//     forAll(curNormalTraction, faceI)
//     {
//         if
//         (
//             (magSqr(this->valueFraction()[faceI]) > 1-SMALL)
//          && (curNormalTraction[faceI] >= law().sigmaMax().value())
//         )
//         {
//             facesToBreakList.insert(faceI);
//             facesToBreakTractionList.insert(curNormalTraction[faceI]);
//         }
//     }

    labelList facesToBreak(facesToBreakList);
    List<scalar> facesToBreakTraction(facesToBreakTractionList);

    label nFacesToBreak = facesToBreak.size();

//     bool breakOnlyOneFace = true;

    label gNFacesToBreak = returnReduce(nFacesToBreak, sumOp<label>());

    if (breakOnlyOneFace_ && gNFacesToBreak)
    {
        if (nFacesToBreak > 1)
        {
            nFacesToBreak = 1;
        }

        // Select face with maximum normal traction
        label faceToBreakIndex = -1;
        scalar faceToBreakTraction = 0;
        forAll(facesToBreakTraction, faceI)
        {
            if (facesToBreakTraction[faceI] > faceToBreakTraction)
            {
                faceToBreakTraction = facesToBreakTraction[faceI];
                faceToBreakIndex = facesToBreak[faceI];
            }
        }

        scalar gMaxTraction =
            returnReduce(faceToBreakTraction, maxOp<scalar>());

        if (Pstream::parRun())
        {
            bool procHasFaceToBreak = false;
            if (nFacesToBreak > 0)
            {
                if ( mag(gMaxTraction - faceToBreakTraction) < SMALL )
                {
                    // Maximum traction is on this processor
                    procHasFaceToBreak = true;
                }
                else
                {
                    nFacesToBreak = 0;
                }
            }

            // Check if maximum is present on more then one processors

            label procID = Pstream::nProcs();
            if (procHasFaceToBreak)
            {
                procID = Pstream::myProcNo();
            }

            label minProcID =
                returnReduce<label>(procID, minOp<label>());

            if (procID != minProcID)
            {
                nFacesToBreak = 0;
            }
        }

        facesToBreak.setSize(nFacesToBreak);
        if (nFacesToBreak)
        {
            facesToBreak[0] = faceToBreakIndex;
        }
    }

    forAll(facesToBreak, fI)
    {
        label faceI = facesToBreak[fI];

        // Switch to full traction boundary condition
        this->valueFraction()[faceI] = symmTensor::zero;
        crazeIndicator_[faceI] = 1;
        crackIndicator_[faceI] = 0;

        if (beta_ > SMALL)
        {
        }
        else if (beta_ < -SMALL)
        {
            initiationTraction_[faceI] =
                n[faceI]*curNormalTraction[faceI]
              + curTangentialTraction[faceI]*law().sigmaMax().value()
               /mag(curTangentialTraction[faceI]);

            traction_[faceI] = initiationTraction_[faceI];
        }
        else
        {
            initiationTraction_[faceI] =
                n[faceI]*law().sigmaMax().value()
              + curTangentialTraction[faceI];

//         initiationTraction_[faceI] =
//             n[faceI]*law().sigmaMax().value()
//           + curTangentialTraction[faceI];

//         // Cohesive traction
//         vector cohesiveTractionIncrement =
//             initiationTraction_[faceI]
//           - oldTraction[faceI];

        Pout << "Crack has started at face: " << faceI
            << ", normal traction: " << curNormalTraction[faceI]
            << ", tangential traction: "
            << curTangentialTraction[faceI]
            << ", initiation traction: "
            << initiationTraction_[faceI]
            << endl;
        }
    }

    this->refGrad() =
        tractionBoundaryGradient().snGrad
        (
            traction_,
            scalarField(patch().size(), 0),
            DName, // field name
            DName, // total field name
            patch(),
            false // not incremental
        );


    gNFacesToBreak = returnReduce(nFacesToBreak, sumOp<label>());

    return gNFacesToBreak;
}


void cohesiveZoneFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (timeSeries_.valid())
    {
        vectorField disp
        (
            patch().size(),
            timeSeries_()(this->db().time().timeOutputValue())
        );

        this->refValue() = disp;
    }

    // Patch normal
    vectorField n = this->patch().nf();

    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        forAll(unloadingSeparationDistance_, faceI)
        {
            vector curSepDist =
                separationDistance_[faceI];

            if (explicitSeparationDistance_)
            {
                curSepDist = oldSeparationDistance_[faceI];
            }

            scalar curNormalSepDist = (n[faceI]&curSepDist);

            if (curNormalSepDist < 0)
            {
                curSepDist -= n[faceI]*(n[faceI]&curSepDist);
                curNormalSepDist = 0;
            }

            vector curTangentialSepDist =
                ((I - n[faceI]*n[faceI])&curSepDist);

            scalar minUnloadingSeparationDistance_ = 0;

            if (beta_ > SMALL) // Mixed mode
            {
                scalar curEffSepDist =
                    sqrt
                    (
                        sqr(curNormalSepDist)
                      + sqr(beta_)*magSqr(curTangentialSepDist)
                    );

                if
                (
                    curEffSepDist
                  > minUnloadingSeparationDistance_*law().deltaC().value()
                )
                {
                    if
                    (
                        curEffSepDist
                      > mag(unloadingSeparationDistance_[faceI])
                    )
                    {
                        unloadingSeparationDistance_[faceI] =
                            curNormalSepDist*n[faceI]
                          + beta_*curTangentialSepDist;
                    }
                }
            }
            else if (beta_ < -SMALL) // Mode II
            {
                scalar curEffSepDist = mag(curTangentialSepDist);

                if
                (
                    curEffSepDist
                  > minUnloadingSeparationDistance_*law().deltaC().value()
                )
                {
                    if
                    (
                        curEffSepDist
                      > mag(unloadingSeparationDistance_[faceI])
                    )
                    {
                        unloadingSeparationDistance_[faceI] =
                            curTangentialSepDist;
                    }
                }
            }
            else // Mode I
            {
                if
                (
                    curNormalSepDist
                  > minUnloadingSeparationDistance_*law().deltaC().value()
                )
                {
                    if
                    (
                        curNormalSepDist
                      > mag(unloadingSeparationDistance_[faceI])
                    )
                    {
                        unloadingSeparationDistance_[faceI] =
                            curNormalSepDist*n[faceI];
                    }
                }
            }


//             if (curSepDist < 0)
//             {
//                 curSepDist = 0;
//             }

//             if
//             (
//                 // Unloading only if
//                 curSepDist > law().deltaC().value()
//             )
//             {
//                 if
//                 (
//                     curSepDist > unloadingSeparationDistance_[faceI]
//                 )
//                 {
//                     unloadingSeparationDistance_[faceI] = curSepDist;
//                 }
//             }
        }

        oldSeparationDistance_ = separationDistance_;

        totRefValue_ = this->refValue();
//         totRefValue_ += this->refValue();

        curTimeIndex_ = this->db().time().timeIndex();
    }

    word DName = this->dimensionedInternalField().name();

    const fvsPatchField<symmTensor>& sigma =
        patch().lookupPatchField<surfaceSymmTensorField, symmTensor>
        (
            "sigmaf"
        );

    // Update componentMixed boundary condition on point displacement field
    pointVectorField& pointD =
        const_cast<pointVectorField&>
        (
            this->db().lookupObject<pointVectorField>("point" + DName)
        );

    if
    (
        pointD.boundaryField()[this->patch().index()].type()
     == componentMixedPointPatchVectorField::typeName
    )
    {
        componentMixedPointPatchVectorField& patchPointD =
            refCast<componentMixedPointPatchVectorField>
            (
                pointD.boundaryField()[this->patch().index()]
            );

        const labelListList& pointFaces = this->patch().patch().pointFaces();
        const vectorField& pointNormals = this->patch().patch().pointNormals();

        patchPointD.valueFraction() = vector::zero;

        if (timeSeries_.valid())
        {
            vectorField disp
            (
                patchPointD.refValue().size(),
                timeSeries_()(this->db().time().timeOutputValue())
            );

            patchPointD.refValue() = disp;
        }

        forAll(pointFaces, pointI)
        {
            const labelList& curPointFaces = pointFaces[pointI];

            bool fixesValue = false;
            forAll(curPointFaces, faceI)
            {
                if (magSqr(valueFraction()[curPointFaces[faceI]]) > SMALL)
                {
                    fixesValue = true;
                    break;
                }

//                 if (magSqr(valueFraction()[curPointFaces[faceI]]) > SMALL)
//                 {
//                     if (patch().type() == wallFvPatch::typeName)
//                     {
//                         patchPointD.valueFraction()[pointI] = vector::one;
//                     }
//                     else
//                     {
//                         patchPointD.valueFraction()[pointI] =
//                             cmptMultiply
//                             (
//                                 pointNormals[pointI],
//                                 pointNormals[pointI]
//                             );
//                     }

//                     break;
//                 }
            }

            if (fixesValue)
            {
                if (patch().type() == wallFvPatch::typeName)
                {
                    patchPointD.valueFraction()[pointI] = vector::one;
                }
                else
                {
                    patchPointD.valueFraction()[pointI] =
                        cmptMultiply
                        (
                            pointNormals[pointI],
                            pointNormals[pointI]
                        );
                }
            }
        }
    }

    // Patch displacement
    const vectorField& D = *this;

    // Patch stress
//     symmTensorField curSigma = sigma;
//     symmTensorField curSigma = sigma + DSigma;

//     // Old traction
//     vectorField oldTraction = (n&sigma);

    // Current traction
    traction_ = (n & sigma);

    // Current normal traction
    scalarField curNormalTraction = (n & (n & sigma));

    // Current tangential traction
    vectorField curTangentialTraction = ((I-n*n)&(n&sigma));

//     // Current normal traction increment
//     scalarField curNormalTractionIncrement = (n&(n&DSigma));

    // Separation distance
    vectorField newSeparationDistance = (totRefValue_ - D);

    if (explicitSeparationDistance_)
    {
        separationDistance_ = newSeparationDistance;
    }
    else
    {

        separationDistance_ =
            separationDistance_
          + relaxationFactor_*(newSeparationDistance - separationDistance_);
    }

    if (patch().type() != wallFvPatch::typeName)
    {
        separationDistance_ *= 2;
    }

    label nCrackedFaces = 0;

    vectorField cohesiveTraction = traction_;
    scalarField cohesivePressure(patch().size(), 0);

    // Chech crack propagation
    forAll(curNormalTraction, faceI)
    {
        vector curSepDist = separationDistance_[faceI];
        if (explicitSeparationDistance_)
        {
            curSepDist = oldSeparationDistance_[faceI];
        }

        scalar curNormalSepDist = (n[faceI]&curSepDist);
//         scalar curNegNormalSepDist = 0;

        if (curNormalSepDist < 0)
        {
            curSepDist -= n[faceI]*(n[faceI]&curSepDist);
//             curNegNormalSepDist = curNormalSepDist;
            curNormalSepDist = 0;
        }

        vector curTangentialSepDist =
            ((I - n[faceI]*n[faceI]) & curSepDist);


        if ((curTangentialSepDist & curTangentialTraction[faceI]) < 0)
        {
            curTangentialSepDist = vector::zero;
        }

//         scalar curNormalTraction = 0;
//         vector curTangentialTraction = vector::zero;


//         if (curSepDist < 0)
//         {
//             curSepDist = 0;
//         }

        if(magSqr(this->valueFraction()[faceI]) < SMALL)
        {
            nCrackedFaces++;

            if (beta_ > SMALL) //Mixed mode
            {
            }
            else if (beta_ < -SMALL) // Mode II
            {
                if
                (
                    (
                        mag(curTangentialSepDist)
                      > law().deltaC().value()
                    )
                 || (
                        mag(curTangentialSepDist)
                     <= mag(unloadingSeparationDistance_[faceI])
                    )
                )
                {
                    // Traction free
                    cohesiveTraction[faceI] = vector::zero;

                    crazeIndicator_[faceI] = 0;
                    crackIndicator_[faceI] = 1;
                }
                else
                {
                    cohesiveTraction =
                        initiationTraction_[faceI]
                       *law().traction(mag(curTangentialSepDist))
                       /law().sigmaMax().value();

                    crazeIndicator_[faceI] = 1;
                    crackIndicator_[faceI] = 0;
                }
            }
            else // Mode I
            {
                if
                (
                    (curNormalSepDist > law().deltaC().value())
                 || (
                        curNormalSepDist
                      < (mag(unloadingSeparationDistance_[faceI])-SMALL)
                    )
                )
                {
                    // Traction free
                    cohesiveTraction[faceI] = vector::zero;

                    crazeIndicator_[faceI] = 0;
                    crackIndicator_[faceI] = 1;

//                 Info << "Traction free face " << faceI
//                     << ", curSepDist: " << curSepDist
//                     << ", deltaC: " << law().deltaC().value()
//                     << ", unloadingSepDist: "
//                     << unloadingSeparationDistance_[faceI] << endl;
                }
                else
                {
                    // Calculate cohesive traction from cohesive zone model

                    cohesiveTraction[faceI] =
                        initiationTraction_[faceI]
                       *law().traction(curNormalSepDist)
                       /law().sigmaMax().value();

                    crazeIndicator_[faceI] = 1;
                    crackIndicator_[faceI] = 0;
                }
            }
        }
    }

    this->refGrad() =
        tractionBoundaryGradient().snGrad
        (
            cohesiveTraction,
            cohesivePressure,
            DName, // field name
            DName, // total field name
            patch(),
            false // not incremental
        );


//     reduce(nCrackedFaces, sumOp<label>());

//     Info << "nCrackFaces: " << nCrackedFaces << endl;

    directionMixedDisplacementFvPatchVectorField::updateCoeffs();
}

// void cohesiveZoneFvPatchVectorField::evaluate
// (
//     const Pstream::commsTypes
// )
// {
//     if (!this->updated())
//     {
//         this->updateCoeffs();
//     }

//     Field<vector> normalValue = transform(valueFraction(), refValue());

//     if (secondOrderCorr_)
//     {
//         word DUName = this->dimensionedInternalField().name();

//         const fvPatchField<tensor>& gradU =
//             patch().lookupPatchField<volTensorField, tensor>
//             (
//                 "grad(" + DUName + ")"
//             );

//         const vectorField n = patch().nf();

//         vectorField nGradUp = (n&gradU.patchInternalField());

//         Field<vector> gradValue =
//             this->patchInternalField()
//           + 0.5*nGradUp/this->patch().deltaCoeffs()
//           + 0.5*refGrad()/this->patch().deltaCoeffs();

//         Field<vector> transformGradValue =
//             transform(I - valueFraction(), gradValue);

//         Field<vector>::operator=(normalValue + transformGradValue);
//     }
//     else
//     {
//         Field<vector> gradValue =
//             this->patchInternalField()
//           + refGrad()/this->patch().deltaCoeffs();

//         Field<vector> transformGradValue =
//             transform(I - valueFraction(), gradValue);

//         Field<vector>::operator=(normalValue + transformGradValue);
//     }

//     transformFvPatchField<vector>::evaluate();
// }


// Write
void cohesiveZoneFvPatchVectorField::write(Ostream& os) const
{
    directionMixedDisplacementFvPatchVectorField::write(os);

    totRefValue_.writeEntry("totRefValue", os);
//     os.writeKeyword("totalFieldName") << totalFieldName_
//         << token::END_STATEMENT << nl;
    os.writeKeyword("simpleCohesiveLaw") << law().type()
        << token::END_STATEMENT << nl;
    crazeIndicator_.writeEntry("crazeIndicator", os);
    crackIndicator_.writeEntry("crackIndicator", os);
    os.writeKeyword("relaxationFactor") << relaxationFactor_
        << token::END_STATEMENT << nl;
    law().writeDict(os);

    os.writeKeyword("explicitSeparationDistance")
        << explicitSeparationDistance_ << token::END_STATEMENT << nl;

    os.writeKeyword("breakOnlyOneFace")
        << breakOnlyOneFace_ << token::END_STATEMENT << nl;

    traction_.writeEntry("traction", os);

    initiationTraction_.writeEntry("initiationTraction", os);

    if (timeSeries_.valid())
    {
        timeSeries_().write(os);
    }

//     os.writeKeyword("nonLinear")
//         << nonLinear_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    cohesiveZoneFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
