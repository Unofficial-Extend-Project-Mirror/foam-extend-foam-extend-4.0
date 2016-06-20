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

\*---------------------------------------------------------------------------*/

#include "cohesiveLawFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
// #include "rheologyModel.H"
#include "regionSplit.H"
#include "crackerFvMesh.H"

#include "solidSolver.H"
// #include "unsULLSSolid.H"
#include "unsTotalLagrangianSolid.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

cohesiveLawFvPatchVectorField::cohesiveLawFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    directionMixedDisplacementFvPatchVectorField(p, iF),
    cohesiveLawPtr_(NULL),
    relaxationFactor_(1.0),
    traction_(p.size(), vector::zero),
    initiationTraction_(p.size(), vector::zero),
    separationDistance_(p.size(), vector::zero),
    oldSeparationDistance_(p.size(), vector::zero),
    unloadingSeparationDistance_(p.size(), vector::zero),
    minUnloadingSeparationDistance_(0.0),
    beta_(0.0),
    explicitSeparationDistance_(false),
    contact_(false),
    contactConstant_(0),
    curTimeIndex_(-1)
{}


cohesiveLawFvPatchVectorField::cohesiveLawFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    directionMixedDisplacementFvPatchVectorField(p, iF),
    cohesiveLawPtr_
    (
        simpleCohesiveLaw::New(dict.lookup("cohesiveLaw"), dict).ptr()
    ),
    relaxationFactor_(readScalar(dict.lookup("relaxationFactor"))),
    traction_(p.size(), vector::zero),
    initiationTraction_(p.size(), vector::zero),
    separationDistance_(p.size(), vector::zero),
    oldSeparationDistance_(p.size(), vector::zero),
    unloadingSeparationDistance_(p.size(), vector::zero),
    minUnloadingSeparationDistance_
    (
        readScalar(dict.lookup("minUnloadingSeparationDistance"))
    ),
    beta_(readScalar(dict.lookup("beta"))),
    explicitSeparationDistance_(dict.lookup("explicitSeparationDistance")),
    contact_(dict.lookup("contact")),
    contactConstant_(readScalar(dict.lookup("contactConstant"))),
    curTimeIndex_(-1)
{
    if (minUnloadingSeparationDistance_ < 0.001)
    {
        minUnloadingSeparationDistance_ = 0.001;
    }

    if (minUnloadingSeparationDistance_ > 1.0)
    {
        minUnloadingSeparationDistance_ = 1.0;
    }

    Info << "minUnloadingSeparationDistance: "
        << minUnloadingSeparationDistance_ << endl;

//     if (beta_ < (1 - SMALL))
//     {
//         beta_ = 0;
//     }
//     else
//     {
//         beta_ = 1;
//     }

    Info << "beta: " << beta_ << endl;

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
        this->valueFraction() = symmTensor::zero;
    }

    if (dict.found("value"))
    {
        vectorField::operator=(vectorField("value", dict, p.size()));
    }
    else
    {
        vectorField normalValue = transform(valueFraction(), refValue());

        vectorField gradValue =
            this->patchInternalField() + refGrad()/this->patch().deltaCoeffs();

        vectorField transformGradValue =
            transform(I - valueFraction(), gradValue);

        vectorField::operator=(normalValue + transformGradValue);
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

    if (dict.found("separationDistance"))
    {
        separationDistance_ =
            vectorField("separationDistance", dict, p.size());
    }

    if (dict.found("oldSeparationDistance"))
    {
        separationDistance_ =
            vectorField("oldSeparationDistance", dict, p.size());
    }

    if (dict.found("unloadingSeparationDistance"))
    {
        unloadingSeparationDistance_ =
            vectorField("unloadingSeparationDistance", dict, p.size());
    }
}


cohesiveLawFvPatchVectorField::cohesiveLawFvPatchVectorField
(
    const cohesiveLawFvPatchVectorField& cpf
)
:
    directionMixedDisplacementFvPatchVectorField(cpf),
    cohesiveLawPtr_(cpf.cohesiveLawPtr_->clone().ptr()),
    relaxationFactor_(cpf.relaxationFactor_),
    traction_(cpf.traction_),
    initiationTraction_(cpf.initiationTraction_),
    separationDistance_(cpf.separationDistance_),
    oldSeparationDistance_(cpf.oldSeparationDistance_),
    unloadingSeparationDistance_(cpf.unloadingSeparationDistance_),
    minUnloadingSeparationDistance_(cpf.minUnloadingSeparationDistance_),
    beta_(cpf.beta_),
    explicitSeparationDistance_(cpf.explicitSeparationDistance_),
    contact_(cpf.contact_),
    contactConstant_(cpf.contactConstant_),
    curTimeIndex_(-1)
{}


cohesiveLawFvPatchVectorField::cohesiveLawFvPatchVectorField
(
    const cohesiveLawFvPatchVectorField& cpf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    directionMixedDisplacementFvPatchVectorField(cpf, p, iF, mapper),
    cohesiveLawPtr_(cpf.cohesiveLawPtr_->clone().ptr()),
    relaxationFactor_(cpf.relaxationFactor_),
    traction_(cpf.traction_, mapper),
    initiationTraction_(cpf.initiationTraction_, mapper),
    separationDistance_(cpf.separationDistance_, mapper),
    oldSeparationDistance_(cpf.oldSeparationDistance_, mapper),
    unloadingSeparationDistance_(cpf.unloadingSeparationDistance_, mapper),
    minUnloadingSeparationDistance_(cpf.minUnloadingSeparationDistance_),
    beta_(cpf.beta_),
    explicitSeparationDistance_(cpf.explicitSeparationDistance_),
    contact_(cpf.contact_),
    contactConstant_(cpf.contactConstant_),
    curTimeIndex_(-1)
{}


cohesiveLawFvPatchVectorField::cohesiveLawFvPatchVectorField
(
    const cohesiveLawFvPatchVectorField& cpf,
    const DimensionedField<vector, volMesh>& iF
)
:
    directionMixedDisplacementFvPatchVectorField(cpf, iF),
    cohesiveLawPtr_(cpf.cohesiveLawPtr_->clone().ptr()),
    relaxationFactor_(cpf.relaxationFactor_),
    traction_(cpf.traction_),
    initiationTraction_(cpf.initiationTraction_),
    separationDistance_(cpf.separationDistance_),
    oldSeparationDistance_(cpf.oldSeparationDistance_),
    unloadingSeparationDistance_(cpf.unloadingSeparationDistance_),
    minUnloadingSeparationDistance_(cpf.minUnloadingSeparationDistance_),
    beta_(cpf.beta_),
    explicitSeparationDistance_(cpf.explicitSeparationDistance_),
    contact_(cpf.contact_),
    contactConstant_(cpf.contactConstant_),
    curTimeIndex_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const simpleCohesiveLaw& cohesiveLawFvPatchVectorField::law() const
{
    if (!cohesiveLawPtr_)
    {
        FatalErrorIn
        (
            "const cohesiveLaw& cohesiveLawFvPatchVectorField::law() const"
        )   << "Law pointer not set" << abort(FatalError);
    }

    return *cohesiveLawPtr_;
}


tmp<scalarField> cohesiveLawFvPatchVectorField::normalTraction() const
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


tmp<scalarField>
cohesiveLawFvPatchVectorField::relativeSeparationDistance() const
{
    tmp<scalarField> tRelativeSeparationDistance
    (
        new scalarField(size(), 0.0)
    );

    vectorField n = patch().nf();

    for (label i=0; i<size(); i++)
    {
        vector curSepDist =
            separationDistance_[i];

        if (mag(n[i]&curSepDist) < 0)
        {
            curSepDist -= n[i]*(n[i]&curSepDist);
        }

        tRelativeSeparationDistance()[i] =
            mag(curSepDist)/law().deltaC().value();
    }

    return tRelativeSeparationDistance;
}


tmp<scalarField> cohesiveLawFvPatchVectorField::effectiveTraction
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
    else
    {
        effTraction() = normalTraction;
    }

    return effTraction;
}


scalar cohesiveLawFvPatchVectorField::effectiveTraction
(
    const vector& traction,
    const vector& normal
) const
{
    scalar normalTraction = (normal&traction);
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
    else
    {
        effTraction = normalTraction;
    }

    return effTraction;
}


void cohesiveLawFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    if (cohesiveLawPtr_ == NULL)
    {
        FatalErrorIn("cohesiveFvPatchVectorField::autoMap")
            << "NULL cohesive law"
                << abort(FatalError);
    }

    directionMixedFvPatchVectorField::autoMap(m);

    traction_.autoMap(m);
    initiationTraction_.autoMap(m);
    separationDistance_.autoMap(m);
    oldSeparationDistance_.autoMap(m);
    unloadingSeparationDistance_.autoMap(m);

    // Must use primitive mesh data because fvMesh data is packed in fields
    // and cannot be accessed during mapping.  HJ, 12/Dec/2008
    vectorField n = patch().patch().faceNormals();

    const labelList& addressing = m.directAddressing();

    label nNewFaces = m.size() - m.sizeBeforeMapping();

    if ( (patch().size()==1) && (nNewFaces == 1) )
    {
        label i=0;
        this->valueFraction()[i] = symmTensor::zero;
        traction_[i] = law().sigmaMax().value()*n[i];
        initiationTraction_[i] = law().sigmaMax().value()*n[i];
        separationDistance_[i] = vector::zero;
        oldSeparationDistance_[i] = vector::zero;
        unloadingSeparationDistance_[i] = vector::zero;
    }
    else if ( (patch().size()==2) && (nNewFaces == 1) )
    {
        label i=1;
        this->valueFraction()[i] = symmTensor::zero;
        traction_[i] = law().sigmaMax().value()*n[i];
        initiationTraction_[i] = law().sigmaMax().value()*n[i];
        separationDistance_[i] = vector::zero;
        oldSeparationDistance_[i] = vector::zero;
        unloadingSeparationDistance_[i] = vector::zero;
    }
    else if ( (patch().size()==2) && (nNewFaces == 2) )
    {
        label i=0;
        this->valueFraction()[i] = symmTensor::zero;
        traction_[i] = law().sigmaMax().value()*n[i];
        initiationTraction_[i] = law().sigmaMax().value()*n[i];
        separationDistance_[i] = vector::zero;
        oldSeparationDistance_[i] = vector::zero;
        unloadingSeparationDistance_[i] = vector::zero;
        i=1;
        this->valueFraction()[i] = symmTensor::zero;
        traction_[i] = law().sigmaMax().value()*n[i];
        initiationTraction_[i] = law().sigmaMax().value()*n[i];
        separationDistance_[i] = vector::zero;
        oldSeparationDistance_[i] = vector::zero;
        unloadingSeparationDistance_[i] = vector::zero;
    }
    else
    {
        for (label i = 1; i < patch().size(); i++)
        {
            if (addressing[i] == 0)
            {
                this->valueFraction()[i] = symmTensor::zero;
                traction_[i] = law().sigmaMax().value()*n[i];
                initiationTraction_[i] = law().sigmaMax().value()*n[i];
                separationDistance_[i] = vector::zero;
                oldSeparationDistance_[i] = vector::zero;
                unloadingSeparationDistance_[i] = vector::zero;
            }
        }
    }

//     label sizeByTwo = patch().size()/2;
//     for (label i = 0; i < sizeByTwo; i++)
//     {
//         if (addressing[i] == addressing[sizeByTwo + i])
//         {
//             this->valueFraction()[i] = symmTensor::zero;
//             this->valueFraction()[sizeByTwo + i] = symmTensor::zero;

//             traction_[i] = law().sigmaMax().value()*n[i];
//             traction_[sizeByTwo + i] =
//                 law().sigmaMax().value()*n[sizeByTwo + i];

//             initiationTraction_[i] = law().sigmaMax().value()*n[i];
//             initiationTraction_[sizeByTwo + i] =
//                 law().sigmaMax().value()*n[sizeByTwo + i];

//             separationDistance_[i] = vector::zero;
//             separationDistance_[sizeByTwo + i] = vector::zero;

//             oldSeparationDistance_[i] = vector::zero;
//             oldSeparationDistance_[sizeByTwo + i] = vector::zero;

//             unloadingSeparationDistance_[i] = vector::zero;
//             unloadingSeparationDistance_[sizeByTwo + i] = vector::zero;
//         }
//     }
}


// Reverse-map the given fvPatchField onto this fvPatchField
void cohesiveLawFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    directionMixedFvPatchVectorField::rmap(ptf, addr);

    const cohesiveLawFvPatchVectorField& dmptf =
        refCast<const cohesiveLawFvPatchVectorField>(ptf);

    // No need to grab the cohesive zone pointer more than once
    if (!cohesiveLawPtr_)
    {
        cohesiveLawPtr_ = dmptf.cohesiveLawPtr_->clone().ptr();
    }

    relaxationFactor_ = dmptf.relaxationFactor_;
    traction_.rmap(dmptf.traction_, addr);
    initiationTraction_.rmap(dmptf.initiationTraction_, addr);
    separationDistance_.rmap(dmptf.separationDistance_, addr);
    oldSeparationDistance_.rmap(dmptf.oldSeparationDistance_, addr);
    unloadingSeparationDistance_.rmap
    (
        dmptf.unloadingSeparationDistance_,
        addr
    );
    minUnloadingSeparationDistance_ = dmptf.minUnloadingSeparationDistance_;
    beta_ = dmptf.beta_;
    contact_ = dmptf.contact_;
    contactConstant_ = dmptf.contactConstant_;
}


// Update the coefficients associated with the patch field
void cohesiveLawFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

//     Info << "updateCoeffs - start" << endl;

    vectorField n = patch().nf();

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

            if (beta_ < SMALL) // Mode I
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
            else // Mixed mode I/II
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
        }

        oldSeparationDistance_ = separationDistance_;
        curTimeIndex_ = this->db().time().timeIndex();
    }

    // Get face cells regions
    const unallocLabelList& faceCells = patch().faceCells();

    const fvMesh& mesh = patch().boundaryMesh().mesh();

    if (!isA<crackerFvMesh>(mesh))
    {
        FatalErrorIn
        (
            "void cohesiveLawFvPatchVectorField::updateCoeffs() const"
        )   << "Mesh should be of type: " << crackerFvMesh::typeName
            << abort(FatalError);
    }

    const crackerFvMesh& crackerMesh = refCast<const crackerFvMesh>(mesh);

    const regionSplit& regions = crackerMesh.regions();

    labelField faceCellRegion(size(), -1);
    forAll(faceCellRegion, faceI)
    {
        faceCellRegion[faceI] = regions[faceCells[faceI]];
    }
    labelField globalFaceCellRegion =
        crackerMesh.globalCrackField(faceCellRegion);

//     // Looking up rheology
//     const rheologyModel& rheology =
//         this->db().objectRegistry::lookupObject<rheologyModel>(rheologyName_);

//     const scalarField mu =
//         rheology.mu()().boundaryField()[patch().index()];
//     const scalarField lambda =
//         rheology.lambda()().boundaryField()[patch().index()];

//     const fvPatchField<tensor>& gradU =
//         patch().lookupPatchField<volTensorField, tensor>("grad(U)");


    // Patch total displacement
    vectorField patchDisplacement = *this;
//         const vectorField& UPatch = *this;

    // // Looking up solid solver
    // const solidSolver& solid =
    //     this->db().objectRegistry::lookupObject<solidSolver>
    //     (
    //         "solidProperties"
    //     );

    // if (isA<solidSolvers::unsULLSSolid>(solid))
    // {
    //     const fvPatchField<vector>& D =
    //         patch().lookupPatchField<volVectorField, vector>("D");

    //     patchDisplacement += D;
    // }

//     // Patch displacement
//     const vectorField& UPatch = *this;





    // Global displacement
    vectorField globalUPatch = crackerMesh.globalCrackField(patchDisplacement);

    // Update separation distance
    vectorField newSeparationDistance(patch().size(), vector::zero);
    const labelList& gcfa = crackerMesh.globalCrackFaceAddressing();
    label globalIndex = crackerMesh.localCrackStart();
    for(label i = 0; i < patch().size(); i++)
    {
        newSeparationDistance[i] =
            globalUPatch[gcfa[globalIndex]]
          - globalUPatch[globalIndex];

        globalIndex++;
    }

//     label sizeByTwo = patch().size()/2;
//     for(label i = 0; i < sizeByTwo; i++)
//     {
//         newSeparationDistance[i] = UPatch[sizeByTwo + i] - UPatch[i];
//         newSeparationDistance[sizeByTwo + i] = -newSeparationDistance[i];
//     }

//     newSeparationDistance =
//         beta_*((I-n*n)&newSeparationDistance)
//       + (n&newSeparationDistance)*n;

    if (explicitSeparationDistance_)
    {
        separationDistance_ = newSeparationDistance;
    }
    else
    {
        separationDistance_ =
            separationDistance_
          + relaxationFactor_
           *(newSeparationDistance - separationDistance_);
    }

//     for(label i=0; i<sizeByTwo; i++)
    globalIndex = crackerMesh.localCrackStart();
    for(label i = 0; i < patch().size(); i++)
    {
        vector curSepDist = separationDistance_[i];

        if (explicitSeparationDistance_)
        {
            curSepDist = oldSeparationDistance_[i];
        }

        scalar curNormalSepDist = (n[i]&curSepDist);
        scalar curNegNormalSepDist = 0;

        if (curNormalSepDist < 0)
        {
            curSepDist -= n[i]*(n[i]&curSepDist);
            curNegNormalSepDist = curNormalSepDist;
            curNormalSepDist = 0;
        }

        vector curTangentialSepDist =
            ((I - n[i]*n[i])&curSepDist);

        scalar curNormalTraction = 0;
        vector curTangentialTraction = vector::zero;

//         if (faceCellRegion[i] == faceCellRegion[sizeByTwo+i])
        if
        (
            globalFaceCellRegion[globalIndex]
         == globalFaceCellRegion[gcfa[globalIndex]]
        )
        {
            if (beta_ > SMALL)
            {
                // Mode II or Mixed mode I/II

                scalar curEffSepDist =
                    sqrt
                    (
                        sqr(curNormalSepDist)
                      + sqr(beta_)*mag(curTangentialSepDist)
                    );

                if
                (
                    curEffSepDist
                  > (mag(unloadingSeparationDistance_[i]) - SMALL)
                )
                {
                    // Loading

                    curNormalTraction = (n[i]&initiationTraction_[i]);
                    curTangentialTraction =
                        ((I - n[i]*n[i])&initiationTraction_[i]);

                    if (curNormalTraction >= 0)
                    {
                        // Tension

                        curNormalTraction *=
                            law().traction(curEffSepDist)
                           /law().sigmaMax().value();

                        curTangentialTraction *=
                            law().traction(curEffSepDist)
                           /law().sigmaMax().value();
                    }
                    else
                    {
                        // Compression

                        curTangentialTraction *=
                            law().traction(curEffSepDist)
                           /law().sigmaMax().value();
                    }
                }
                else
                {
                    // Unloading

                    scalar unloadingNormalTraction =
                        (n[i]&initiationTraction_[i]);
                    vector unloadingTangentialTraction =
                        ((I - n[i]*n[i])&initiationTraction_[i]);

                    if (unloadingNormalTraction >= 0)
                    {
                        unloadingNormalTraction *=
                            law().traction
                            (
                                mag(unloadingSeparationDistance_[i])
                            )
                           /law().sigmaMax().value();

                        unloadingTangentialTraction *=
                            law().traction
                            (
                                mag(unloadingSeparationDistance_[i])
                            )
                           /law().sigmaMax().value();

                        curNormalTraction =
                            unloadingNormalTraction
                           *(
                                curEffSepDist
                               /mag(unloadingSeparationDistance_[i])
                            );

                        curTangentialTraction =
                            unloadingTangentialTraction
                           *(
                                curEffSepDist
                               /mag(unloadingSeparationDistance_[i])
                            );
                    }
                    else
                    {
                        unloadingTangentialTraction *=
                            law().traction
                            (
                                mag(unloadingSeparationDistance_[i])
                            )
                           /law().sigmaMax().value();

                        curTangentialTraction =
                            unloadingTangentialTraction
                           *(
                                curEffSepDist
                               /mag(unloadingSeparationDistance_[i])
                            );

                        curNormalTraction = unloadingNormalTraction;
                    }
                }
            }
            else
            {
                // Mode I
                if
                (
                    curNormalSepDist
                  > (mag(unloadingSeparationDistance_[i]) - SMALL)
                )
                {
                    //Loading
                    curNormalTraction =
                        (n[i]&initiationTraction_[i]);

                    curNormalTraction *=
                        law().traction(curNormalSepDist)
                       /law().sigmaMax().value();

                    curTangentialTraction =
                        ((I - n[i]*n[i])&initiationTraction_[i]);

                    curTangentialTraction *=
                        law().traction(curNormalSepDist)
                       /law().sigmaMax().value();
                }
                else
                {
                    // Unloading

                    scalar unloadingNormalTraction =
                        (n[i]&initiationTraction_[i])
                       *law().traction(mag(unloadingSeparationDistance_[i]))
                       /law().sigmaMax().value();

                    vector unloadingTangentialTraction =
                        ((I - n[i]*n[i])&initiationTraction_[i])
                       *law().traction(mag(unloadingSeparationDistance_[i]))
                       /law().sigmaMax().value();

                    curNormalTraction =
                        unloadingNormalTraction
                       *(
                           curNormalSepDist
                          /mag(unloadingSeparationDistance_[i])
                        );

//                     if (curNegNormalSepDist < 0)
//                     {
//                         curNormalTraction =
//                             unloadingNormalTraction
//                            *(
//                                 curNegNormalSepDist
//                                /mag(unloadingSeparationDistance_[i])
//                             );
//                     }

                    curTangentialTraction =
                        unloadingTangentialTraction
                       *(
                           curNormalSepDist
                          /mag(unloadingSeparationDistance_[i])
                        );
                }
            }
        }

        // Correct direction of tangential traction
        // to be aligned with tangential separation distance
        if (mag(curTangentialSepDist) > SMALL)
        {
            curTangentialTraction =
                mag(curTangentialTraction)
               *curTangentialSepDist
               /mag(curTangentialSepDist);
        }

        // Contact
        if(contact_ && (curNegNormalSepDist < 0))
        {
            scalar m =
                (1.0/(1.0-contactConstant_))
               *law().sigmaMax().value()
               /law().deltaC().value();

            curNormalTraction = m*curNegNormalSepDist;

            Info << "Contact: " << curNegNormalSepDist/law().deltaC().value()
                << ", " << curNormalTraction << endl;
        }

        traction_[i] = curNormalTraction*n[i] + curTangentialTraction;
//         traction_[sizeByTwo + i] = -traction_[i];
    }

//     Info << this->refGrad() << endl;
//     Info << traction_.size() << ", " << mu.size() << ", "
//         << lambda.size() << ", " << n.size() << ", "
//         << gradDf.size() << endl;

    word DName = this->dimensionedInternalField().name();

    const fvsPatchField<tensor>& gradDf =
        patch().lookupPatchField<surfaceTensorField, tensor>
        (
            "grad" + DName + "f"
        );

    const fvsPatchField<scalar>& mu =
        patch().lookupPatchField<surfaceScalarField, scalar>
        (
            "muf"
        );

    const fvsPatchField<scalar>& lambda =
        patch().lookupPatchField<surfaceScalarField, scalar>
        (
            "lambdaf"
        );

    this->refGrad() =
    (
        traction_
      - (n & (mu*gradDf.T() - (mu + lambda)*gradDf))
      - n*lambda*tr(gradDf)
    )/(2.0*mu + lambda);

    directionMixedFvPatchVectorField::updateCoeffs();

//     Info << "updateCoeffs - end" << endl;
}


// Write
void cohesiveLawFvPatchVectorField::write(Ostream& os) const
{
    directionMixedFvPatchVectorField::write(os);
    traction_.writeEntry("traction", os);
    initiationTraction_.writeEntry("initiationTraction", os);
    separationDistance_.writeEntry("separationDistance", os);
    oldSeparationDistance_.writeEntry("oldSeparationDistance", os);
    unloadingSeparationDistance_.writeEntry("unloadingSeparationDistance", os);
    os.writeKeyword("cohesiveLaw") << law().type()
        << token::END_STATEMENT << nl;
    os.writeKeyword("relaxationFactor") << relaxationFactor_
        << token::END_STATEMENT << nl;
    os.writeKeyword("beta") << beta_
        << token::END_STATEMENT << nl;
    os.writeKeyword("explicitSeparationDistance")
        << explicitSeparationDistance_
            << token::END_STATEMENT << nl;
    os.writeKeyword("contact") << contact_
        << token::END_STATEMENT << nl;
    os.writeKeyword("contactConstant") << contactConstant_
        << token::END_STATEMENT << nl;
    os.writeKeyword("minUnloadingSeparationDistance")
        << minUnloadingSeparationDistance_
            << token::END_STATEMENT << nl;
    law().writeDict(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, cohesiveLawFvPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
