/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
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

Class
    cohesiveZoneIncrementalFvPatchVectorField

Description

\*---------------------------------------------------------------------------*/

#include "cohesiveZoneIncrementalFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "transformField.H"
#include "volFields.H"
#include "pointFields.H"
#include "constitutiveModel.H"
#include "solidSolver.H"
#include "wallFvPatch.H"
#include "componentMixedPointPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

cohesiveZoneIncrementalFvPatchVectorField
::cohesiveZoneIncrementalFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    directionMixedDisplacementFvPatchVectorField(p, iF),
    totRefValue_(p.size(), vector::zero),
    totalFieldName_("undefined"),
    cohesiveLawPtr_(NULL),
    crackIndicator_(p.size(), 0.0),
    crazeIndicator_(p.size(), 0.0),
    relaxationFactor_(1.0),
    separationDistance_(p.size(), 0.0),
    oldSeparationDistance_(p.size(), 0.0),
    unloadingSeparationDistance_(p.size(), 0.0),
    explicitSeparationDistance_(false),
    curTimeIndex_(-1),
    initiationTraction_(p.size(), vector::zero),
    breakOnlyOneFace_(false)
{}


cohesiveZoneIncrementalFvPatchVectorField
::cohesiveZoneIncrementalFvPatchVectorField
(
    const cohesiveZoneIncrementalFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    directionMixedDisplacementFvPatchVectorField(ptf, p, iF, mapper),
    totRefValue_(ptf.totRefValue_),
    totalFieldName_(ptf.totalFieldName_),
    cohesiveLawPtr_(ptf.cohesiveLawPtr_),
    crackIndicator_(ptf.crackIndicator_),
    crazeIndicator_(ptf.crazeIndicator_),
    relaxationFactor_(ptf.relaxationFactor_),
    separationDistance_(ptf.separationDistance_),
    oldSeparationDistance_(ptf.oldSeparationDistance_),
    unloadingSeparationDistance_(ptf.unloadingSeparationDistance_),
    explicitSeparationDistance_(ptf.explicitSeparationDistance_),
    curTimeIndex_(ptf.curTimeIndex_),
    initiationTraction_(ptf.initiationTraction_),
    breakOnlyOneFace_(ptf.breakOnlyOneFace_)
{}


cohesiveZoneIncrementalFvPatchVectorField
::cohesiveZoneIncrementalFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    directionMixedDisplacementFvPatchVectorField(p, iF),
    totRefValue_(p.size(), vector::zero),
    totalFieldName_(dict.lookup("totalFieldName")),
    cohesiveLawPtr_
    (
        simpleCohesiveLaw::New(dict.lookup("simpleCohesiveLaw"), dict).ptr()
    ),
    crackIndicator_(p.size(), 0.0),
    crazeIndicator_(p.size(), 0.0),
    relaxationFactor_(readScalar(dict.lookup("relaxationFactor"))),
    separationDistance_(p.size(), 0.0),
    oldSeparationDistance_(p.size(), 0.0),
    unloadingSeparationDistance_(p.size(), 0.0),
    explicitSeparationDistance_(dict.lookup("explicitSeparationDistance")),
    curTimeIndex_(-1),
    initiationTraction_(p.size(), vector::zero),
    breakOnlyOneFace_(dict.lookup("breakOnlyOneFace"))
{
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
            scalarField("separationDistance", dict, p.size());
    }

    if (dict.found("oldSeparationDistance"))
    {
        separationDistance_ =
            scalarField("oldSeparationDistance", dict, p.size());
    }
}


cohesiveZoneIncrementalFvPatchVectorField
::cohesiveZoneIncrementalFvPatchVectorField
(
    const cohesiveZoneIncrementalFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    directionMixedDisplacementFvPatchVectorField(ptf, iF),
    totRefValue_(ptf.totRefValue_),
    totalFieldName_(ptf.totalFieldName_),
    cohesiveLawPtr_(ptf.cohesiveLawPtr_),
    crackIndicator_(ptf.crackIndicator_),
    crazeIndicator_(ptf.crazeIndicator_),
    relaxationFactor_(ptf.relaxationFactor_),
    separationDistance_(ptf.separationDistance_),
    oldSeparationDistance_(ptf.oldSeparationDistance_),
    unloadingSeparationDistance_(ptf.unloadingSeparationDistance_),
    explicitSeparationDistance_(ptf.explicitSeparationDistance_),
    curTimeIndex_(ptf.curTimeIndex_),
    initiationTraction_(ptf.initiationTraction_),
    breakOnlyOneFace_(ptf.breakOnlyOneFace_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Map from self
void cohesiveZoneIncrementalFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    if (cohesiveLawPtr_ == NULL)
    {
        FatalErrorIn("cohesiveZoneIncrementalFvPatchVectorField::autoMap")
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
    initiationTraction_.autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void cohesiveZoneIncrementalFvPatchVectorField::rmap
(
    const fvPatchField<vector>& ptf,
    const labelList& addr
)
{
    directionMixedDisplacementFvPatchVectorField::rmap(ptf, addr);

    const cohesiveZoneIncrementalFvPatchVectorField& dmptf =
        refCast<const cohesiveZoneIncrementalFvPatchVectorField>(ptf);

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
    initiationTraction_.rmap(dmptf.initiationTraction_, addr);
    breakOnlyOneFace_ = dmptf.breakOnlyOneFace_;
}


label cohesiveZoneIncrementalFvPatchVectorField::updateCrack()
{
    // Looking up solid solver
    const solidSolver& stress =
        this->db().objectRegistry::lookupObject<solidSolver>
        (
            "solidProperties"
        );

    Switch nonLinear
    (
        stress.solidProperties().lookup("nonLinear")
    );

    Switch enforceLinear
    (
        stress.solidProperties().lookup("enforceLinear")
    );

    // Looking up rheology
    const fvPatchField<scalar>& mu =
        patch().lookupPatchField<volScalarField, scalar>("mu");

    const fvPatchField<scalar>& lambda =
        patch().lookupPatchField<volScalarField, scalar>("lambda");

    word DDName = this->dimensionedInternalField().name();

    const fvPatchField<tensor>& gradDD =
        patch().lookupPatchField<volTensorField, tensor>
        (
            "grad(" + DDName + ")"
        );

    const fvPatchField<tensor>& gradD =
        patch().lookupPatchField<volTensorField, tensor>
        (
            "grad(" + totalFieldName_ + ")"
        );

    const fvPatchField<symmTensor>& sigma =
        patch().lookupPatchField<volSymmTensorField, symmTensor>
        (
            "sigma"
        );

    const fvPatchField<symmTensor>& DSigma =
        patch().lookupPatchField<volSymmTensorField, symmTensor>
        (
            "DSigma"
        );

    // Patch stress
    symmTensorField curSigma = sigma + DSigma;

    // Patch normal
    vectorField n = this->patch().nf();

    // Old traction
    vectorField oldTraction = (n&sigma);

    // Current normal traction
    scalarField curNormalTraction = (n&(n&curSigma));

    // Current tangential traction
    vectorField curTangentialTraction = ((I-n*n)&(n&curSigma));

    // Current normal traction increment
    scalarField curNormalTractionIncrement = (n&(n&DSigma));

    // Select potential faces to breake
    SLList<label> facesToBreakList;
    SLList<scalar> facesToBreakTractionList;

    forAll(curNormalTraction, faceI)
    {
        if
        (
            (magSqr(this->valueFraction()[faceI]) > 1-SMALL)
         && (curNormalTraction[faceI] >= law().sigmaMax().value())
        )
        {
            facesToBreakList.insert(faceI);
            facesToBreakTractionList.insert(curNormalTraction[faceI]);
        }
    }

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

        // Select internal face with maximum normal traction
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

        initiationTraction_[faceI] =
            n[faceI]*law().sigmaMax().value();

//         initiationTraction_[faceI] =
//             n[faceI]*law().sigmaMax().value()
//           + curTangentialTraction[faceI];

        // Cohesive traction
        vector cohesiveTractionIncrement =
            initiationTraction_[faceI]
          - oldTraction[faceI];

        Pout << "Crack has started at face: " << faceI
            << ", normal traction: " << curNormalTraction[faceI]
            << ", tangential traction: "
            << curTangentialTraction[faceI]
            << ", initiation traction: "
            << initiationTraction_[faceI]
            << endl;

        this->refGrad()[faceI] =
        (
            cohesiveTractionIncrement
          - (
                n[faceI]
              & (
                    mu[faceI]*gradDD[faceI].T()
                  - (mu[faceI] + lambda[faceI])*gradDD[faceI]
                )
            )
          - n[faceI]*lambda[faceI]*tr(gradDD[faceI])
        );

        if(nonLinear && !enforceLinear)
        {
            this->refGrad()[faceI] -=
                (
                    n[faceI]
                  & (mu[faceI]*(gradDD[faceI] & gradDD[faceI].T()))
                )
              + (
                    n[faceI]
                  & (mu[faceI]*(gradDD[faceI] & gradD[faceI].T()))
                )
              + (
                    n[faceI]
                  & (mu[faceI]*(gradD[faceI] & gradDD[faceI].T()))
                )
              + 0.5*n[faceI]*lambda[faceI]
               *tr(gradDD[faceI] & gradDD[faceI].T())
              + 0.5*n[faceI]*lambda[faceI]
               *tr(gradDD[faceI] & gradD[faceI].T())
              + 0.5*n[faceI]*lambda[faceI]
               *tr(gradD[faceI] & gradDD[faceI].T());
        }

        if (stress.rheology().plasticityActive())
        {
            this->refGrad()[faceI] +=
                2*mu[faceI]
               *(
                    n[faceI]
                  & stress.rheology()
                   .DEpsilonP().boundaryField()[patch().index()][faceI]
                );
        }

        this->refGrad()[faceI] /=
            (2.0*mu[faceI] + lambda[faceI]);
    }

    gNFacesToBreak = returnReduce(nFacesToBreak, sumOp<label>());

    return gNFacesToBreak;
}


void cohesiveZoneIncrementalFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        forAll(unloadingSeparationDistance_, faceI)
        {
            scalar curSepDist =
                separationDistance_[faceI];

            if (explicitSeparationDistance_)
            {
                curSepDist = oldSeparationDistance_[faceI];
            }

            if (curSepDist < 0)
            {
                curSepDist = 0;
            }

            if
            (
                // Unloading only if
                curSepDist > law().deltaC().value()
            )
            {
                if
                (
                    curSepDist > unloadingSeparationDistance_[faceI]
                )
                {
                    unloadingSeparationDistance_[faceI] = curSepDist;
                }
            }
        }

        oldSeparationDistance_ = separationDistance_;

        totRefValue_ += this->refValue();

        curTimeIndex_ = this->db().time().timeIndex();
    }

    // Looking up solid solver
    const solidSolver& stress =
        this->db().objectRegistry::lookupObject<solidSolver>
        (
            "solidProperties"
        );

    Switch nonLinear
    (
        stress.solidProperties().lookup("nonLinear")
    );

    Switch enforceLinear
    (
        stress.solidProperties().lookup("enforceLinear")
    );

    // Looking up rheology
    const fvPatchField<scalar>& mu =
        patch().lookupPatchField<volScalarField, scalar>("mu");

    const fvPatchField<scalar>& lambda =
        patch().lookupPatchField<volScalarField, scalar>("lambda");

    word DDName = this->dimensionedInternalField().name();

    const fvPatchField<tensor>& gradDD =
        patch().lookupPatchField<volTensorField, tensor>
        (
            "grad(" + DDName + ")"
        );

    const fvPatchField<tensor>& gradD =
        patch().lookupPatchField<volTensorField, tensor>
        (
            "grad(" + totalFieldName_ + ")"
        );

    const fvPatchField<vector>& D =
        patch().lookupPatchField<volVectorField, vector>
        (
            totalFieldName_
        );

    const fvPatchField<symmTensor>& sigma =
        patch().lookupPatchField<volSymmTensorField, symmTensor>
        (
            "sigma"
        );

    const fvPatchField<symmTensor>& DSigma =
        patch().lookupPatchField<volSymmTensorField, symmTensor>
        (
            "DSigma"
        );

    // Update componentMixed boundary condition on point displacement field
    pointVectorField& pointDD =
        const_cast<pointVectorField&>
        (
            this->db().lookupObject<pointVectorField>("point" + DDName)
        );

    if
    (
        pointDD.boundaryField()[this->patch().index()].type()
     == componentMixedPointPatchVectorField::typeName
    )
    {
        componentMixedPointPatchVectorField& patchPointDD =
            refCast<componentMixedPointPatchVectorField>
            (
                pointDD.boundaryField()[this->patch().index()]
            );

        const labelListList& pointFaces = this->patch().patch().pointFaces();
        const vectorField& pointNormals = this->patch().patch().pointNormals();

        patchPointDD.valueFraction() = vector::zero;

        forAll(pointFaces, pointI)
        {
            const labelList& curPointFaces = pointFaces[pointI];

            forAll(curPointFaces, faceI)
            {
                if (magSqr(valueFraction()[curPointFaces[faceI]]) > SMALL)
                {
                    if (patch().type() == wallFvPatch::typeName)
                    {
                        patchPointDD.valueFraction()[pointI] = vector::one;
                    }
                    else
                    {
                        patchPointDD.valueFraction()[pointI] =
                            cmptMultiply
                            (
                                pointNormals[pointI],
                                pointNormals[pointI]
                            );
                    }

                    break;
                }
            }
        }
    }

    // Patch displacement increment
    const vectorField& DD = *this;

    // Patch stress
    symmTensorField curSigma = sigma + DSigma;

    // Patch normal
    vectorField n = this->patch().nf();

    // Old traction
    vectorField oldTraction = (n&sigma);

    // Current normal traction
    scalarField curNormalTraction = (n&(n&curSigma));

    // Current normal traction increment
    scalarField curNormalTractionIncrement = (n&(n&DSigma));

    // Separation distance
    if (explicitSeparationDistance_)
    {
        separationDistance_ = (n&(totRefValue_ - (D + DD)));
    }
    else
    {
        scalarField newSeparationDistance = (n&(totRefValue_ - (D + DD)));

        separationDistance_ =
            separationDistance_
          + relaxationFactor_*(newSeparationDistance - separationDistance_);
    }

    if (patch().type() != wallFvPatch::typeName)
    {
        separationDistance_ *= 2;
    }

//     label nCrackedFaces = 0;

    // Chech crack propagation
    forAll(curNormalTraction, faceI)
    {
        vector cohesiveTractionIncrement = vector::zero;

        scalar curSepDist = separationDistance_[faceI];
        if (explicitSeparationDistance_)
        {
            curSepDist = oldSeparationDistance_[faceI];
        }

        if (curSepDist < 0)
        {
            curSepDist = 0;
        }

        if(magSqr(this->valueFraction()[faceI]) < SMALL)
        {
//             nCrackedFaces++;

            if
            (
                (curSepDist > law().deltaC().value())
             || (curSepDist < (unloadingSeparationDistance_[faceI]-SMALL))
            )
            {
                // Traction free
                cohesiveTractionIncrement =
                    vector::zero
                  - oldTraction[faceI];

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
//                 cohesiveTractionIncrement =
//                     law().traction(curSepDist)*n[faceI]
//                   - oldRelaxedTraction[faceI];

                cohesiveTractionIncrement =
                    initiationTraction_[faceI]
                   *law().traction(curSepDist)/law().sigmaMax().value()
                  - oldTraction[faceI];

                crazeIndicator_[faceI] = 1;
                crackIndicator_[faceI] = 0;
            }

            this->refGrad()[faceI] =
            (
                cohesiveTractionIncrement
              - (
                    n[faceI]
                  & (
                        mu[faceI]*gradDD[faceI].T()
                      - (mu[faceI] + lambda[faceI])*gradDD[faceI]
                    )
                )
              - n[faceI]*lambda[faceI]*tr(gradDD[faceI])
            );

            if(nonLinear && !enforceLinear)
            {
                this->refGrad()[faceI] -=
                    (
                        n[faceI]
                      & (mu[faceI]*(gradDD[faceI] & gradDD[faceI].T()))
                    )
                  + (n[faceI] & (mu[faceI]*(gradDD[faceI] & gradD[faceI].T())))
                  + (n[faceI] & (mu[faceI]*(gradD[faceI] & gradDD[faceI].T())))
                  + 0.5*n[faceI]*lambda[faceI]
                   *tr(gradDD[faceI] & gradDD[faceI].T())
                  + 0.5*n[faceI]*lambda[faceI]
                   *tr(gradDD[faceI] & gradD[faceI].T())
                  + 0.5*n[faceI]*lambda[faceI]
                   *tr(gradD[faceI] & gradDD[faceI].T());
            }

            if (stress.rheology().plasticityActive())
            {
                this->refGrad()[faceI] +=
                    2*mu[faceI]
                   *(
                        n[faceI]
                      & stress.rheology()
                       .DEpsilonP().boundaryField()[patch().index()][faceI]
                    );
            }

            this->refGrad()[faceI] /= (2.0*mu[faceI] + lambda[faceI]);
        }
    }

//     reduce(nCrackedFaces, sumOp<label>());

//     Info << "nCrackFaces: " << nCrackedFaces << endl;

    directionMixedDisplacementFvPatchVectorField::updateCoeffs();
}

// void cohesiveZoneIncrementalFvPatchVectorField::evaluate
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
void cohesiveZoneIncrementalFvPatchVectorField::write(Ostream& os) const
{
    directionMixedDisplacementFvPatchVectorField::write(os);

    totRefValue_.writeEntry("totRefValue", os);
    os.writeKeyword("totalFieldName") << totalFieldName_
        << token::END_STATEMENT << nl;
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

    initiationTraction_.writeEntry("initiationTraction", os);

//     os.writeKeyword("nonLinear")
//         << nonLinear_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    cohesiveZoneIncrementalFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
