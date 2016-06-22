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

#include "solidCohesiveFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "constitutiveModel.H"
#include "regionSplit.H"
#include "crackerFvMesh.H"
#include "tractionBoundaryGradient.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

solidCohesiveFvPatchVectorField::solidCohesiveFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    directionMixedFvPatchVectorField(p, iF),
    fieldName_("undefined"),
    relaxationFactor_(1.0),
    traction_(p.size(), vector::zero),
    //minUnloadingSeparationDistance_(0.0),
    contact_(false),
    cracked_(0,false),
    curTractionN_(0, 0.0),
    oldTractionN_(0, 0.0),
    curTractionS_(0, 0.0),
    oldTractionS_(0, 0.0),
    deltaN_(0, 0.0),
    oldDeltaN_(0, 0.0),
    deltaS_(0, 0.0),
    oldDeltaS_(0, 0.0),
    unloadingDeltaEff_(0, 0.0),
    currentGI_(0, 0.0),
    oldGI_(0, 0.0),
    currentGII_(0, 0.0),
    oldGII_(0, 0.0),
    curTimeIndex_(-1),
    penaltyFactorPtr_(NULL),
    penaltyScale_(0.0),
    frictionCoeff_(0.0),
    explicitSeparationDistance_(false),
    curDeltaN_(0, 0.0),
    curDeltaS_(0, 0.0),
    updateGlobalPatchMaterials_(true),
    globalPatchMaterials_(0, 0.0),
    nonLinear_(nonLinearGeometry::OFF),
    orthotropic_(false)
{}


solidCohesiveFvPatchVectorField::solidCohesiveFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    directionMixedFvPatchVectorField(p, iF),
    fieldName_(dimensionedInternalField().name()),
    relaxationFactor_(readScalar(dict.lookup("relaxationFactor"))),
    traction_(p.size(), vector::zero),
    contact_(dict.lookup("contact")),
    cracked_(p.size(),false),
    curTractionN_(p.size(), 0.0),
    oldTractionN_(p.size(), 0.0),
    curTractionS_(p.size(), 0.0),
    oldTractionS_(p.size(), 0.0),
    deltaN_(p.size(), 0.0),
    oldDeltaN_(p.size(), 0.0),
    deltaS_(p.size(), 0.0),
    oldDeltaS_(p.size(), 0.0),
    unloadingDeltaEff_(p.size(), 0.0),
    currentGI_(p.size(), 0.0),
    oldGI_(p.size(), 0.0),
    currentGII_(p.size(), 0.0),
    oldGII_(p.size(), 0.0),
    curTimeIndex_(-1),
    penaltyFactorPtr_(NULL),
    penaltyScale_(readScalar(dict.lookup("penaltyScale"))),
    frictionCoeff_(readScalar(dict.lookup("frictionCoeff"))),
    explicitSeparationDistance_(dict.lookup("explicitSeparationDistance")),
    curDeltaN_(p.size(), 0.0),
    curDeltaS_(p.size(), 0.0),
    updateGlobalPatchMaterials_(true),
    globalPatchMaterials_(0, 0.0),
    nonLinear_(nonLinearGeometry::OFF),
    orthotropic_(false)
{
    Info<< "Creating solidCohesive patch" << nl
        << "\tOnly Dugdale law currently available!" << endl;

    // Check if traction boundary is for non linear solver
    if (dict.found("nonLinear"))
    {
        nonLinear_ = nonLinearGeometry::nonLinearNames_.read
        (
            dict.lookup("nonLinear")
        );

        if (nonLinear_ == nonLinearGeometry::UPDATED_LAGRANGIAN)
        {
            Info << "\tnonLinear set to updated Lagrangian"
                << endl;
        }
        else if (nonLinear_ == nonLinearGeometry::TOTAL_LAGRANGIAN)
        {
            Info << "\tnonLinear set to total Lagrangian"
                << endl;
        }
    }

    if (dict.found("orthotropic"))
    {
        orthotropic_ = Switch(dict.lookup("orthotropic"));
        Info << "\t\torthotropic set to " << orthotropic_ << endl;
    }

    // if (minUnloadingSeparationDistance_ < 0.001)
    // {
    //     minUnloadingSeparationDistance_ = 0.001;
    // }

    // if (minUnloadingSeparationDistance_ > 1.0)
    // {
    //     minUnloadingSeparationDistance_ = 1.0;
    // }

    // Info << "minUnloadingSeparationDistance: "
    //     << minUnloadingSeparationDistance_ << endl;

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

    // if (dict.found("unloadingSeparationDistance"))
    // {
    //     unloadingSeparationDistance_ =
    //         vectorField("unloadingSeparationDistance", dict, p.size());
    // }

    if (dict.found("cracked"))
    {
//      cracked_ = Field<bool>("cracked", dict, p.size());
        cracked_ = Field<scalar>("cracked", dict, p.size());
    }

    if (dict.found("curTractionN"))
    {
        curTractionN_ = scalarField("curTractionN", dict, p.size());
    }
    if (dict.found("curTractionS"))
    {
        curTractionS_ = scalarField("curTractionS", dict, p.size());
    }
    if (dict.found("oldTractionN"))
    {
        oldTractionN_ = scalarField("oldTractionN", dict, p.size());
    }
    if (dict.found("oldTractionS"))
    {
        oldTractionS_ = scalarField("oldTractionS", dict, p.size());
    }

    if (dict.found("deltaN"))
    {
        deltaN_ = scalarField("deltaN", dict, p.size());
    }
    if (dict.found("deltaS"))
    {
        deltaS_ = scalarField("deltaS", dict, p.size());
    }
    if (dict.found("oldDeltaN"))
    {
        oldDeltaN_ = scalarField("oldDeltaN", dict, p.size());
    }
    if (dict.found("oldDeltaS"))
    {
        oldDeltaS_ = scalarField("oldDeltaS", dict, p.size());
    }
    if (dict.found("curDeltaN"))
    {
        curDeltaN_ = scalarField("curDeltaN", dict, p.size());
    }
    if (dict.found("curDeltaS"))
    {
        curDeltaS_ = scalarField("curDeltaS", dict, p.size());
    }
    if (dict.found("unloadingDeltaEff"))
    {
        unloadingDeltaEff_ = scalarField("unloadingDeltaEff", dict, p.size());
    }

    if (dict.found("currentGI"))
    {
        currentGI_ = scalarField("currentGI", dict, p.size());
    }
    if (dict.found("currentGII"))
    {
        currentGII_ = scalarField("currentGII", dict, p.size());
    }
    if (dict.found("oldGI"))
    {
        oldGI_ = scalarField("oldGI", dict, p.size());
    }
    if (dict.found("oldGII"))
    {
        oldGII_ = scalarField("oldGII", dict, p.size());
    }
}


solidCohesiveFvPatchVectorField::solidCohesiveFvPatchVectorField
(
    const solidCohesiveFvPatchVectorField& cpf
)
:
    directionMixedFvPatchVectorField(cpf),
    fieldName_(cpf.fieldName_),
    relaxationFactor_(cpf.relaxationFactor_),
    traction_(cpf.traction_),
    //minUnloadingSeparationDistance_(cpf.minUnloadingSeparationDistance_),
    contact_(cpf.contact_),
    cracked_(cpf.cracked_),
    curTractionN_(cpf.curTractionN_),
    oldTractionN_(cpf.oldTractionN_),
    curTractionS_(cpf.curTractionS_),
    oldTractionS_(cpf.oldTractionS_),
    deltaN_(cpf.deltaN_),
    oldDeltaN_(cpf.oldDeltaN_),
    deltaS_(cpf.deltaS_),
    oldDeltaS_(cpf.oldDeltaS_),
    unloadingDeltaEff_(cpf.unloadingDeltaEff_),
    currentGI_(cpf.currentGI_),
    oldGI_(cpf.oldGI_),
    currentGII_(cpf.currentGII_),
    oldGII_(cpf.oldGII_),
    curTimeIndex_(-1),
    penaltyFactorPtr_(cpf.penaltyFactorPtr_),
    penaltyScale_(cpf.penaltyScale_),
    frictionCoeff_(cpf.frictionCoeff_),
    explicitSeparationDistance_(cpf.explicitSeparationDistance_),
    curDeltaN_(cpf.curDeltaN_),
    curDeltaS_(cpf.curDeltaS_),
    updateGlobalPatchMaterials_(cpf.updateGlobalPatchMaterials_),
    globalPatchMaterials_(cpf.globalPatchMaterials_),
    nonLinear_(cpf.nonLinear_),
    orthotropic_(cpf.orthotropic_)
{}


solidCohesiveFvPatchVectorField::solidCohesiveFvPatchVectorField
(
    const solidCohesiveFvPatchVectorField& cpf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    directionMixedFvPatchVectorField(cpf, p, iF, mapper),
    fieldName_(cpf.fieldName_),
    relaxationFactor_(cpf.relaxationFactor_),
    traction_(cpf.traction_, mapper),
    //minUnloadingSeparationDistance_(cpf.minUnloadingSeparationDistance_),
    contact_(cpf.contact_),
    cracked_(cpf.cracked_),
    curTractionN_(cpf.curTractionN_),
    oldTractionN_(cpf.oldTractionN_),
    curTractionS_(cpf.curTractionS_),
    oldTractionS_(cpf.oldTractionS_),
    deltaN_(cpf.deltaN_),
    oldDeltaN_(cpf.oldDeltaN_),
    deltaS_(cpf.deltaS_),
    oldDeltaS_(cpf.oldDeltaS_),
    unloadingDeltaEff_(cpf.unloadingDeltaEff_),
    currentGI_(cpf.currentGI_),
    oldGI_(cpf.oldGI_),
    currentGII_(cpf.currentGII_),
    oldGII_(cpf.oldGII_),
    curTimeIndex_(-1),
    penaltyFactorPtr_(cpf.penaltyFactorPtr_),
    penaltyScale_(cpf.penaltyScale_),
    frictionCoeff_(cpf.frictionCoeff_),
    explicitSeparationDistance_(cpf.explicitSeparationDistance_),
    curDeltaN_(cpf.curDeltaN_),
    curDeltaS_(cpf.curDeltaS_),
    updateGlobalPatchMaterials_(cpf.updateGlobalPatchMaterials_),
    globalPatchMaterials_(cpf.globalPatchMaterials_),
    nonLinear_(cpf.nonLinear_),
    orthotropic_(cpf.orthotropic_)
{}


solidCohesiveFvPatchVectorField::solidCohesiveFvPatchVectorField
(
    const solidCohesiveFvPatchVectorField& cpf,
    const DimensionedField<vector, volMesh>& iF
)
:
    directionMixedFvPatchVectorField(cpf, iF),
    fieldName_(cpf.fieldName_),
    relaxationFactor_(cpf.relaxationFactor_),
    traction_(cpf.traction_),
    //minUnloadingSeparationDistance_(cpf.minUnloadingSeparationDistance_),
    contact_(cpf.contact_),
    cracked_(cpf.cracked_),
    curTractionN_(cpf.curTractionN_),
    oldTractionN_(cpf.oldTractionN_),
    curTractionS_(cpf.curTractionS_),
    oldTractionS_(cpf.oldTractionS_),
    deltaN_(cpf.deltaN_),
    oldDeltaN_(cpf.oldDeltaN_),
    deltaS_(cpf.deltaS_),
    oldDeltaS_(cpf.oldDeltaS_),
    unloadingDeltaEff_(cpf.unloadingDeltaEff_),
    currentGI_(cpf.currentGI_),
    oldGI_(cpf.oldGI_),
    currentGII_(cpf.currentGII_),
    oldGII_(cpf.oldGII_),
    curTimeIndex_(-1),
    penaltyFactorPtr_(cpf.penaltyFactorPtr_),
    penaltyScale_(cpf.penaltyScale_),
    frictionCoeff_(cpf.frictionCoeff_),
    explicitSeparationDistance_(cpf.explicitSeparationDistance_),
    curDeltaN_(cpf.curDeltaN_),
    curDeltaS_(cpf.curDeltaS_),
    updateGlobalPatchMaterials_(cpf.updateGlobalPatchMaterials_),
    globalPatchMaterials_(cpf.globalPatchMaterials_),
    nonLinear_(cpf.nonLinear_),
    orthotropic_(cpf.orthotropic_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<scalarField> solidCohesiveFvPatchVectorField::crackingAndDamage() const
{
    tmp<scalarField> tcrackingAndDamage
    (
        new scalarField(size(), 1.0)
    );
    scalarField& cad = tcrackingAndDamage();

    forAll(cad, facei)
    {
        if (cracked_[facei])
        {
            cad[facei] = 2.0;
        }
    }

    return tcrackingAndDamage;
}


tmp<scalarField> solidCohesiveFvPatchVectorField::GI() const
{
    tmp<scalarField> tGI
    (
        new scalarField(size(), 0.0)
    );
    scalarField& GI = tGI();

    forAll(GI, facei)
    {
        GI[facei] = currentGI_[facei];
    }

    return tGI;
}


tmp<scalarField> solidCohesiveFvPatchVectorField::GII() const
{
    tmp<scalarField> tGII
    (
        new scalarField(size(), 0.0)
    );
    scalarField& GII = tGII();

    forAll(GII, facei)
    {
        GII[facei] = currentGII_[facei];
    }

    return tGII;
}


bool solidCohesiveFvPatchVectorField::cracking()
{
    const fvMesh& mesh = patch().boundaryMesh().mesh();
    const crackerFvMesh& crackerMesh = refCast<const crackerFvMesh>(mesh);

    // global fields
    Field<scalar> globalCracked =
        crackerMesh.globalCrackField(cracked_);

    label sumDamaged = 0;
    label sumCracked = 0;

    forAll(globalCracked, facei)
    {
        if (globalCracked[facei] > 0.0)
        {
            sumCracked++;
        }
        else
        {
            sumDamaged++;
        }
    }
    Info<< "\t\tThere are " << sumDamaged << " damaged faces" << nl
        << "\t\tThere are " << sumCracked << " cracked faces" << endl;

    return false;
}


void solidCohesiveFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    directionMixedFvPatchVectorField::autoMap(m);

    traction_.autoMap(m);
    cracked_.autoMap(m);
    curTractionN_.autoMap(m);
    oldTractionN_.autoMap(m);
    curTractionS_.autoMap(m);
    oldTractionS_.autoMap(m);
    deltaN_.autoMap(m);
    oldDeltaN_.autoMap(m);
    deltaS_.autoMap(m);
    oldDeltaS_.autoMap(m);
    unloadingDeltaEff_.autoMap(m);
    currentGI_.autoMap(m);
    oldGI_.autoMap(m);
    currentGII_.autoMap(m);
    oldGII_.autoMap(m);
    curDeltaN_.autoMap(m);
    curDeltaS_.autoMap(m);
    globalPatchMaterials_.autoMap(m);

    // force this field to be reset next time updateCoeffs is called
    updateGlobalPatchMaterials_ = true;

    // Must use primitive mesh data because fvMesh data is packed in fields
    // and cannot be accessed during mapping.  HJ, 12/Dec/2008
    vectorField n = patch().patch().faceNormals();

    const labelList& addressing = m.directAddressing();

    label nNewFaces = m.size() - m.sizeBeforeMapping();

    // philipc:
    // when a new face is created in the crack, we need to set
    // its traction, but how do we do it:
    // we calculate the traciton in the solver (updatedCrack.H)
    // but we cannot access that field here during mapping
    // So in the solver we will have to cast this patch and update
    // the traction_ field after mapping and call updateCoeffs.
    // For here, we will just set traction to zero.

    if ( (patch().size() == 1) && (nNewFaces == 1) )
    {
        label i=0;
        this->valueFraction()[i] = symmTensor::zero;

        traction_[i] = vector::zero;
        cracked_[i] = false;
        curTractionN_[i] = 0.0;
        oldTractionN_[i] = 0.0;
        curTractionS_[i] = 0.0;
        oldTractionS_[i] = 0.0;
        deltaN_[i] = 0.0;
        oldDeltaN_[i] = 0.0;
        deltaS_[i] = 0.0;
        oldDeltaS_[i] = 0.0;
        curDeltaN_[i] = 0.0;
        curDeltaS_[i] = 0.0;
        unloadingDeltaEff_[i] = 0.0;
        currentGI_[i] = 0.0;
        oldGI_[i] = 0.0;
        currentGII_[i] = 0.0;
        oldGII_[i] = 0.0;
    }
    else if ( (patch().size() == 2) && (nNewFaces == 1) )
    {
        label i=1;
        this->valueFraction()[i] = symmTensor::zero;
        traction_[i] = vector::zero;
        cracked_[i] = false;

        curTractionN_[i] = 0.0;
        oldTractionN_[i] = 0.0;
        curTractionS_[i] = 0.0;
        oldTractionS_[i] = 0.0;
        deltaN_[i] = 0.0;
        oldDeltaN_[i] = 0.0;
        deltaS_[i] = 0.0;
        oldDeltaS_[i] = 0.0;
        curDeltaN_[i] = 0.0;
        curDeltaS_[i] = 0.0;
        unloadingDeltaEff_[i] = 0.0;
        currentGI_[i] = 0.0;
        oldGI_[i] = 0.0;
        currentGII_[i] = 0.0;
        oldGII_[i] = 0.0;
    }
    else if ( (patch().size()==2) && (nNewFaces == 2) )
    {
        label i=0;
        this->valueFraction()[i] = symmTensor::zero;
        traction_[i] = vector::zero;
        i=1;
        this->valueFraction()[i] = symmTensor::zero;
        traction_[i] = vector::zero;

        cracked_[i] = false;
        curTractionN_[i] = 0.0;
        oldTractionN_[i] = 0.0;
        curTractionS_[i] = 0.0;
        oldTractionS_[i] = 0.0;
        deltaN_[i] = 0.0;
        oldDeltaN_[i] = 0.0;
        deltaS_[i] = 0.0;
        oldDeltaS_[i] = 0.0;
        curDeltaN_[i] = 0.0;
        curDeltaS_[i] = 0.0;
        unloadingDeltaEff_[i] = 0.0;
        currentGI_[i] = 0.0;
        oldGI_[i] = 0.0;
        currentGII_[i] = 0.0;
        oldGII_[i] = 0.0;
    }
    else
    {
        for (label i = 1; i < patch().size(); i++)
        {
            if (addressing[i] == 0)
            {
                this->valueFraction()[i] = symmTensor::zero;
                traction_[i] = vector::zero;

                cracked_[i] = false;
                curTractionN_[i] = 0.0;
                oldTractionN_[i] = 0.0;
                curTractionS_[i] = 0.0;
                oldTractionS_[i] = 0.0;
                deltaN_[i] = 0.0;
                oldDeltaN_[i] = 0.0;
                deltaS_[i] = 0.0;
                oldDeltaS_[i] = 0.0;
                curDeltaN_[i] = 0.0;
                curDeltaS_[i] = 0.0;
                unloadingDeltaEff_[i] = 0.0;
                currentGI_[i] = 0.0;
                oldGI_[i] = 0.0;
                currentGII_[i] = 0.0;
                oldGII_[i] = 0.0;
            }
        }
    }
}


// Reverse-map the given fvPatchField onto this fvPatchField
void solidCohesiveFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    directionMixedFvPatchVectorField::rmap(ptf, addr);

    const solidCohesiveFvPatchVectorField& dmptf =
        refCast<const solidCohesiveFvPatchVectorField>(ptf);

    relaxationFactor_ = dmptf.relaxationFactor_;
    traction_.rmap(dmptf.traction_, addr);
    contact_ = dmptf.contact_;
    cracked_ = dmptf.cracked_;
    curTractionN_ = dmptf.curTractionN_;
    oldTractionN_ = dmptf.oldTractionN_;
    curTractionS_ = dmptf.curTractionS_;
    oldTractionS_ = dmptf.oldTractionS_;
    deltaN_ = dmptf.deltaN_;
    oldDeltaN_ = dmptf.oldDeltaN_;
    deltaS_ = dmptf.deltaS_;
    oldDeltaS_ = dmptf.oldDeltaS_;
    curDeltaN_ = dmptf.curDeltaN_;
    curDeltaS_ = dmptf.curDeltaS_;
    globalPatchMaterials_ = dmptf.globalPatchMaterials_;
    unloadingDeltaEff_ = dmptf.unloadingDeltaEff_;
    currentGI_ = dmptf.currentGI_;
    currentGII_ = dmptf.currentGII_;
}


// Update the coefficients associated with the patch field
void solidCohesiveFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    vectorField n = patch().nf();

    // lookup cohesive properties from rheology
    const constitutiveModel& rheology =
        this->db().objectRegistry::lookupObject<constitutiveModel>
        (
            "rheologyProperties"
        );

    label patchID = patch().index();
    const scalarField sigmaMax =
        rheology.cohLaw().sigmaMax()().boundaryField()[patchID];
    const scalarField tauMax =
        rheology.cohLaw().tauMax()().boundaryField()[patchID];
    const scalarField GIc =
        rheology.cohLaw().GIc()().boundaryField()[patchID];
    const scalarField GIIc =
        rheology.cohLaw().GIIc()().boundaryField()[patchID];

    // update old values
    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        // we force the penalty factor to be calculated here for the first time
        // as all processors must call this at the same time
        if (contact_ && patch().size() > 0)
        {
            // force calculation of penalty factor here
            penaltyFactor();
        }

        oldDeltaN_ = deltaN_;
        oldDeltaS_ = deltaS_;
        // curDelta is always latest delta
        // whereas delta is only latest when explicitDist is false
        //oldDeltaN_ = curDeltaN_;
        //oldDeltaS_ = curDeltaS_;
        oldTractionN_ = curTractionN_;
        oldTractionS_ = curTractionS_;
        oldGI_ = currentGI_;
        oldGII_ = currentGII_;
        unloadingDeltaEff_ =
            max
            (
                unloadingDeltaEff_,
                Foam::sqrt
                (
                    sqr(max(deltaN_, scalar(0))) + sqr(deltaS_)
                )
            );

        curTimeIndex_ = this->db().time().timeIndex();
    }


    // Get face cells regions
    const unallocLabelList& faceCells = patch().faceCells();
    const fvMesh& mesh = patch().boundaryMesh().mesh();

    if (!isA<crackerFvMesh>(mesh))
    {
        FatalErrorIn
        (
            "void solidCohesiveFvPatchVectorField::updateCoeffs() const"
        )   << "Mesh should be of type: " << crackerFvMesh::typeName
            << abort(FatalError);
    }

    const crackerFvMesh& crackerMesh = refCast<const crackerFvMesh>(mesh);

    // global patch material field is needed for multimaterial cohesive laws
    if (updateGlobalPatchMaterials_)
    {
        updateGlobalPatchMaterials_ = false;

        if (mesh.objectRegistry::foundObject<volScalarField>("materials"))
        {
            scalarField localPatchMaterials =
                patch().lookupPatchField<volScalarField, scalar>
                ("materials").patchInternalField();
            scalarField newGlobalPatchMaterials =
                crackerMesh.globalCrackField(localPatchMaterials);

            globalPatchMaterials_.setSize(newGlobalPatchMaterials.size());
            globalPatchMaterials_ = newGlobalPatchMaterials;
        }
    }

    const regionSplit& regions = crackerMesh.regions();

    labelField faceCellRegion(size(), -1);

    forAll(faceCellRegion, faceI)
    {
        faceCellRegion[faceI] = regions[faceCells[faceI]];
    }

    labelField globalFaceCellRegion =
        crackerMesh.globalCrackField(faceCellRegion);

    // Patch displacement
    vectorField UPatch = *this;
    if (fieldName_ == "DU")
    {
        UPatch += patch().lookupPatchField<volVectorField, vector>("U");
    }

    // Global displacement
    vectorField globalUPatch = crackerMesh.globalCrackField(UPatch);


    // Update separation distance
    vectorField delta(patch().size(), vector::zero);
    const labelList& gcfa = crackerMesh.globalCrackFaceAddressing();
    label globalIndex = crackerMesh.localCrackStart();
    for (label i = 0; i < patch().size(); i++)
    {
        // newSeparationDistance[i] =
        delta[i] =
            globalUPatch[gcfa[globalIndex]]
          - globalUPatch[globalIndex];

        globalIndex++;
    }

    //globalIndex = crackerMesh.localCrackStart();
    for (label i = 0; i < patch().size(); i++)
    {
        // update deltas
        curDeltaN_[i] = n[i] & delta[i];
        curDeltaS_[i] = mag( (I - sqr(n[i])) & delta[i]); // shearing
        if (explicitSeparationDistance_)
        {
            deltaN_[i] = oldDeltaN_[i];
            deltaS_[i] = oldDeltaS_[i];

            if (mag(deltaN_[i]) < SMALL)
            {
                deltaN_[i] = 2*SMALL;
            }
            FatalError
                << "explicitSeparationDistance_ is true!" << abort(FatalError);
        }
        else
        {
            // relax
            deltaN_[i] =
                relaxationFactor_*curDeltaN_[i]
              + (1.0 - relaxationFactor_)*deltaN_[i];

            deltaS_[i] =
                relaxationFactor_*curDeltaS_[i]
              + (1.0 - relaxationFactor_)*deltaS_[i];
        }

        // update tractions
        curTractionN_[i] = n[i] & traction_[i];
        curTractionS_[i] = mag( (I-sqr(n[i])) & traction_[i] );

        // effective delta - only positive deltaN is considered
        const scalar deltaEff =
            Foam::sqrt(max(deltaN_[i], 0.0)*max(deltaN_[i], 0.0)
          + deltaS_[i]*deltaS_[i]);

        // update energies
        // stop calculating after cracking for convergence
        // because crack might jump in and out of damaged/failed
        // only update energies if there is loading
        if ( !cracked_[i] && (deltaEff > (unloadingDeltaEff_[i]-SMALL)) )
        {
            // if the average normal stress is tensile
            if ((curTractionN_[i]+oldTractionN_[i]) > 0.0)
            {
                // trapezoidal rule
                currentGI_[i] =
                    oldGI_[i]
                  + ((0.5*(curTractionN_[i]+oldTractionN_[i]))*
                    (deltaN_[i]-oldDeltaN_[i]));
            }
            else
            {
                // no modeI energy lost if in compression
                currentGI_[i] = oldGI_[i];
            }

            // mode II - trapezoidal rule
            currentGII_[i] = oldGII_[i]
              + ((0.5*(curTractionS_[i]+oldTractionS_[i]))*
                (deltaS_[i]-oldDeltaS_[i]));
        }
        //scalar currentG = currentGI_[i] + currentGII_[i];

        // if
        // (
        //  ( globalFaceCellRegion[globalIndex]
        //    == globalFaceCellRegion[gcfa[globalIndex]] )
        //  ||
        //  !cracked_[i]
        //  )

        // new tractions to be calculated
        scalar curNormalTraction = 0;
        vector curTangentialTraction = vector::zero;

        // loading
        // if the effective delta is greater than unloadingDeltaEff then
        // there is loading
        // unloadingDeltaEff is the maximum previsouly reached deltaEff
        if (deltaEff > (unloadingDeltaEff_[i]-SMALL) )
        {
            // at the moment loading and unloading are the same
            // if total energy is dissipated, then fully crack face
            //if ( currentG > GIc[i] || cracked_[i] ) //Gc )
            // propagation
            if ( ((currentGI_[i]/GIc[i]) + (currentGII_[i]/GIIc[i])) >= 1 )
            {
                //Pout << "GIc[i] is " << GIc[i] << ", curG is " << currentG << endl;
                if (!cracked_[i])
                {
                    Pout << "Face " << i << " is fully cracked" << endl;
                }

                cracked_[i] = true;


                // failed faces might come in to contact so we need to deal with them
                // here we could use the face delta values but that assumes that the
                // opposing crack faces are aligned direclty opposite one another, which
                // in general they are not. So we actually need to calculate the actual
                // face penetration distances and then a standard penalty approach
                // contact is probably fine. We will use simple assumption for now which
                // is fine is the relative displacement
                // of the opposing faces is not too large
                // To-do: we must calculate actual distances
                curNormalTraction = 0.0;
                curTangentialTraction = vector::zero;
                if ( contact_ && deltaN_[i] <= 0.0 )
                {
                    curNormalTraction = deltaN_[i]*penaltyFactor();
                    //Info << "penaltyFactor() is " << penaltyFactor() << endl;

                    // friction
                    vector sDir = (I - sqr(n[i]))&delta[i];
                    sDir /= mag(sDir+vector(SMALL,SMALL,SMALL));
                    scalar slip = mag(deltaS_[gcfa[i]] - deltaS_[i]);
                    scalar fricMag =
                        min
                        (
                            slip*penaltyFactor(),
                            frictionCoeff_*curNormalTraction // stick/slip
                        );
                    curTangentialTraction = fricMag*sDir;
                }

                // relax tractions
                curNormalTraction =
                    relaxationFactor_*curNormalTraction
                  + (1-relaxationFactor_)*(n[i] & traction_[i]);

                curTangentialTraction =
                    relaxationFactor_*curTangentialTraction
                  + (1-relaxationFactor_)*( (I -sqr(n[i])) & traction_[i]);

                traction_[i] = curNormalTraction*n[i] + curTangentialTraction;
            }
            // damging face with positive normal delta
            else if ( deltaN_[i] > 0.0 )
            {
                if (cracked_[i])
                {
                    Pout << "Face " << i << " is un-cracked" << endl;
                }

                cracked_[i] = false;

                // set traction in a fixed point iteration manner to force
                // (tN/sigmaMax)^2 + (tS/tauMax)^2 = 1
                // while also allowing varying mode mix by assuming colinear traction
                // and delta
                // Dugdale version
                // scalar tN =
                //   sigmaMax[i] * deltaN_[i] /
                //   (SMALL + ::sqrt( (deltaN_[i]*deltaN_[i])
                // + (deltaS_[i]*deltaS_[i])*(sigmaMax[i]*sigmaMax[i]
                // /(tauMax[i]*tauMax[i])) ));
                // scalar tS =
                //   tauMax[i] * deltaS_[i] /
                //   (SMALL + ::sqrt( (deltaS_[i]*deltaS_[i])
                // + (deltaN_[i]*deltaN_[i])*(tauMax[i]*tauMax[i]
                // /(sigmaMax[i]*sigmaMax[i])) ));

                // Given current deltaN and deltaS, the cohesive law
                // (Dugdale, linear, etc)gives back current traction:
                scalar tN = n[i] & traction_[i];
                scalar tS = mag( (I - sqr(n[i])) & traction_[i] );
                // Update tN and tS
                rheology.cohLaw().damageTractions
                (
                    tN,
                    tS,
                    deltaN_[i],
                    deltaS_[i],
                    currentGI_[i],
                    currentGII_[i],
                    i,
                    globalPatchMaterials_
                );

                // shear delta direction
                vector sDir = (I - sqr(n[i]))&delta[i];
                sDir /= mag(sDir+vector(SMALL,SMALL,SMALL));

                // relax tractions
                curNormalTraction =
                    relaxationFactor_*tN
                  + (1-relaxationFactor_)*(n[i] & traction_[i]);

                curTangentialTraction =
                    relaxationFactor_*(tS*sDir)
                  + (1-relaxationFactor_)*( (I -sqr(n[i])) & traction_[i]);

                traction_[i] = curNormalTraction*n[i] + curTangentialTraction;
            }
            // damaging faces with negative normal delta
            else
            {
                if (cracked_[i])
                {
                    Pout << "Face " << i << " is un-cracked" << endl;
                }

                cracked_[i] = false;
                //Pout << "Contact and shearing face " << i << endl;

                //Pout << "Face " << i << " is damaging but in contact" << endl;
                // set shear traction from cohesive law
                // set normal traction to contact condition
                scalar tS = tauMax[i];

                vector sDir = (I - sqr(n[i])) & delta[i];
                sDir /= mag(sDir+vector(SMALL,SMALL,SMALL));

                // Simple penalty condition
                scalar penaltyFac = 0.0;
                if (contact_)
                {
                    penaltyFac = penaltyFactor();
                }
                scalar contactTN = deltaN_[i]*penaltyFac;

                // relax tractions
                curNormalTraction =
                    relaxationFactor_*contactTN
                  + (1-relaxationFactor_)*(n[i] & traction_[i]);

                curTangentialTraction =
                    relaxationFactor_*(tS*sDir)
                  + (1-relaxationFactor_)*( (I -sqr(n[i])) & traction_[i]);

                traction_[i] = curNormalTraction*n[i] + curTangentialTraction;
            }
        }
        else
        {
            // unloading
            // as the current effective delta is less than the old time
            // effective delta

            // We have two choice for unloading:
            // (a) ductile approach
            //         unload with initial stiffness (which is infinite)
            //         similar to ductilve metal stress-strain curve
            //         as we use infinite intial stiffness, this means to elastic
            //         energy is recovered
            //         and we immediately start loading in opposite direction
            // (b) brittle approach
            //         unload straight back to the origin
            //         this means we recover elastic energy
            //
            // approach (b) is numerically "nicer"; however, for Dugdale cohesive
            // zone this implys that only half the energy is dissipated just before
            // failure and then the other half is dissipated at failure. This may
            // not be an issue in practive but it does not really make sense.
            // Obviously it is fine for a linear cohesive zone as the energy
            // is smoothly dissipated up to failure. For now, we will implement
            // approach (b), but this requires more thinking...

            // reduce traction linearly with the reduction in delta
            scalar scaleFactor = deltaEff/(unloadingDeltaEff_[i]);
            traction_[i] =
                relaxationFactor_*scaleFactor*traction_[i]
              + (1 - relaxationFactor_)*traction_[i];
        }
    }

    bool incremental(fieldName_ == "DU");

    this->refGrad() = tractionBoundaryGradient::snGrad
    (
        traction_,
        scalarField(traction_.size(), 0.0),
        fieldName_,
        "U",
        patch(),
        orthotropic_,
        nonLinear_,
        incremental
    );

    directionMixedFvPatchVectorField::updateCoeffs();
}


void solidCohesiveFvPatchVectorField::calcPenaltyFactor()
{
    // calculate penalty factor similar to standardPenalty contact model
    // approx penaltyFactor from mechanical properties
    // this can then be scaled using the penaltyScale
    // to-do: write equivalent for orthotropic
    if (orthotropic_)
    {
        FatalErrorIn
        (
            "void solidCohesiveFvPatchVectorField::calcPenaltyFactor()"
        )   << "solidCohesiveFvPatchVectorField::calcPenaltyFactor()"
            << " has yet to be written for orthotropic"
            << exit(FatalError);
    }

    const label patchID = patch().index();
    const fvMesh& mesh = patch().boundaryMesh().mesh();

    const scalarField& mu =
        mesh.objectRegistry::lookupObject<volScalarField>
        ("mu").boundaryField()[patchID];

    const scalarField& lambda =
        mesh.objectRegistry::lookupObject<volScalarField>
        ("lambda").boundaryField()[patchID];

    // avarage contact patch bulk modulus
    scalar bulkModulus = gAverage(lambda + (2.0/3.0)*mu);

    // average contact patch face area
    scalar faceArea = gAverage(mesh.magSf().boundaryField()[patchID]);

    // average contact patch cell volume
    scalar cellVolume = 0.0;

    const volScalarField::DimensionedInternalField & V = mesh.V();
    {
        const unallocLabelList& faceCells =
            mesh.boundary()[patchID].faceCells();

        forAll(mesh.boundary()[patchID], facei)
        {
            cellVolume += V[faceCells[facei]];
        }
    }
    cellVolume /= patch().size();

    // approximate penalty factor based on Hallquist et al.
    // we approximate penalty factor for traction instead of force
    penaltyFactorPtr_ =
        new scalar(penaltyScale_*bulkModulus*faceArea/cellVolume);
}


// Write
void solidCohesiveFvPatchVectorField::write(Ostream& os) const
{
    directionMixedFvPatchVectorField::write(os);
    traction_.writeEntry("traction", os);
    os.writeKeyword("relaxationFactor")
        << relaxationFactor_ << token::END_STATEMENT << nl;

    cracked_.writeEntry("cracked", os);
    curTractionN_.writeEntry("curTractionN", os);
    curTractionS_.writeEntry("curTractionS", os);
    oldTractionN_.writeEntry("oldTractionN", os);
    oldTractionS_.writeEntry("oldTractionS", os);
    deltaN_.writeEntry("deltaN", os);
    deltaS_.writeEntry("deltaS", os);
    oldDeltaN_.writeEntry("oldDeltaN", os);
    oldDeltaS_.writeEntry("oldDeltaS", os);
    unloadingDeltaEff_.writeEntry("unloadingDeltaEff", os);
    currentGI_.writeEntry("currentGI", os);
    currentGII_.writeEntry("currentGII", os);
    oldGI_.writeEntry("oldGI", os);
    oldGII_.writeEntry("oldGII", os);

    os.writeKeyword("curTimeIndex")
        << curTimeIndex_ << token::END_STATEMENT << nl;
    os.writeKeyword("contact") << contact_<< token::END_STATEMENT << nl;
    os.writeKeyword("penaltyScale")
        << penaltyScale_ << token::END_STATEMENT << nl;
    os.writeKeyword("nonLinear")
        << nonLinearGeometry::nonLinearNames_[nonLinear_]
        << token::END_STATEMENT << nl;
    os.writeKeyword("orthotropic")
        << orthotropic_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, solidCohesiveFvPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
