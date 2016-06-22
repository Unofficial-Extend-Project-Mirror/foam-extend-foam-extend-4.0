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
    solidContactFvPatchVectorField

\*---------------------------------------------------------------------------*/

#include "solidContactFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "transformField.H"
#include "volFields.H"
#include "tractionBoundaryGradient.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

solidContactFvPatchVectorField::solidContactFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    directionMixedFvPatchVectorField(p, iF),
    fieldName_("undefined"),
    master_("undefined"),
    contactActive_(false),
    rigidMaster_(false),
    normalContactModelPtr_(NULL),
    frictionContactModelPtr_(NULL),
    shadowPatchID_(-1),
    patchInChargeOfCorrection_(-1),
    masterFaceZoneName_("undefined"),
    slaveFaceZoneName_("undefined"),
    masterFaceZoneID_(-1),
    slaveFaceZoneID_(-1),
    // masterFaceZonePatchPoints_(pointField(0)),
    // slaveFaceZonePatchPoints_(pointField(0)),
    // masterFaceZonePatchFaces_(faceList(0)),
    // slaveFaceZonePatchFaces_(faceList(0)),
    masterFaceZonePatchPtr_(NULL),
    slaveFaceZonePatchPtr_(NULL),
    interpolationMethod_("undefined"),
    slaveToMasterPatchToPatchInterpolatorPtr_(NULL),
    slaveToMasterGgiInterpolatorPtr_(NULL),
    masterFaceZonePatchInterpolatorPtr_(NULL),
    slaveFaceZonePatchInterpolatorPtr_(NULL),
    oldMasterFaceZonePoints_(0,vector::zero),
    oldSlaveFaceZonePoints_(0,vector::zero),
    alg_(Foam::intersection::VISIBLE),
    dir_(Foam::intersection::CONTACT_SPHERE),
    curTimeIndex_(-1),
    iCorr_(0),
    correctionFreq_(0),
    orthotropic_(false),
    stickSlipFieldPtr_(NULL),
    forceCorrection_(false),
    nonLinear_(nonLinearGeometry::OFF)
{}


solidContactFvPatchVectorField::solidContactFvPatchVectorField
(
    const solidContactFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    directionMixedFvPatchVectorField(ptf, p, iF, mapper),
    fieldName_(ptf.fieldName_),
    master_(ptf.master_),
    contactActive_(ptf.contactActive_),
    rigidMaster_(ptf.rigidMaster_),
    normalContactModelPtr_(ptf.normalContactModelPtr_),
    frictionContactModelPtr_(ptf.frictionContactModelPtr_),
    shadowPatchID_(ptf.shadowPatchID_),
    patchInChargeOfCorrection_(ptf.patchInChargeOfCorrection_),
    masterFaceZoneName_(ptf.masterFaceZoneName_),
    slaveFaceZoneName_(ptf.slaveFaceZoneName_),
    masterFaceZoneID_(ptf.masterFaceZoneID_),
    slaveFaceZoneID_(ptf.slaveFaceZoneID_),
    masterFaceZonePatchPtr_(ptf.masterFaceZonePatchPtr_),
    slaveFaceZonePatchPtr_(ptf.slaveFaceZonePatchPtr_),
    interpolationMethod_(ptf.interpolationMethod_),
    slaveToMasterPatchToPatchInterpolatorPtr_
    (
        ptf.slaveToMasterPatchToPatchInterpolatorPtr_
    ),
    slaveToMasterGgiInterpolatorPtr_(ptf.slaveToMasterGgiInterpolatorPtr_),
    masterFaceZonePatchInterpolatorPtr_
    (
        ptf.masterFaceZonePatchInterpolatorPtr_
    ),
    slaveFaceZonePatchInterpolatorPtr_(ptf.slaveFaceZonePatchInterpolatorPtr_),
    oldMasterFaceZonePoints_(ptf.oldMasterFaceZonePoints_),
    oldSlaveFaceZonePoints_(ptf.oldSlaveFaceZonePoints_),
    alg_(ptf.alg_),
    dir_(ptf.dir_),
    curTimeIndex_(ptf.curTimeIndex_),
    iCorr_(ptf.iCorr_),
    correctionFreq_(ptf.correctionFreq_),
    orthotropic_(ptf.orthotropic_),
    stickSlipFieldPtr_(ptf.stickSlipFieldPtr_),
    forceCorrection_(ptf.forceCorrection_),
    nonLinear_(ptf.nonLinear_)
{}


solidContactFvPatchVectorField::solidContactFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
    : // only the master reads the properties
    directionMixedFvPatchVectorField(p, iF),
    fieldName_(dimensionedInternalField().name()),
    master_
    (
        checkPatchAndFaceZones(dict)
        ? dict.lookup("master") : dict.lookup("master")
        ),
    contactActive_(dict.lookup("contactActive")),
    rigidMaster_(false), //dict.lookup("rigidMaster")),
    normalContactModelPtr_(NULL),
    frictionContactModelPtr_(NULL),
    shadowPatchID_
    (
        patch().patch().boundaryMesh().findPatchID(dict.lookup("shadowPatch"))
        ),
    patchInChargeOfCorrection_
    ( bool(patch().index()<shadowPatchID_) ? true : false
        ),
    masterFaceZoneName_
    (
        master_
        ? patch().name() + "FaceZone"
        : patch().boundaryMesh().mesh().boundary()[shadowPatchID_].name()
        + "FaceZone"
        ),
    slaveFaceZoneName_
    (
        master_
        ? patch().boundaryMesh().mesh().boundary()[shadowPatchID_].name()
        + "FaceZone"
        : patch().name() + "FaceZone"
        ),
    masterFaceZoneID_
    (
        patch().boundaryMesh().mesh().faceZones().findZoneID
        (
            masterFaceZoneName_
            )
        ),
    slaveFaceZoneID_
    (
        patch().boundaryMesh().mesh().faceZones().findZoneID(slaveFaceZoneName_)
        ),
    masterFaceZonePatchPtr_
    (
         master_ ?
         new PrimitivePatch<face, Foam::List, pointField>
         (
             patch().boundaryMesh().mesh().faceZones()[masterFaceZoneID_]()
             .localFaces(),
             patch().boundaryMesh().mesh().faceZones()[masterFaceZoneID_]()
                 .localPoints()
      )
     :
     NULL
     ),
     slaveFaceZonePatchPtr_
     (
         master_ ?
     new PrimitivePatch<face, Foam::List, pointField>
         (
             patch().boundaryMesh().mesh().faceZones()[slaveFaceZoneID_]()
             .localFaces(),
             patch().boundaryMesh().mesh().faceZones()[slaveFaceZoneID_]()
             .localPoints()
      )
     :
     NULL
     ),
     interpolationMethod_
     (
         master_ ?
     dict.lookup("interpolationMethod")
     :
     word("interpolationMethodUndefinedForSlave")
     ),
     slaveToMasterPatchToPatchInterpolatorPtr_(NULL),
     slaveToMasterGgiInterpolatorPtr_(NULL),
     masterFaceZonePatchInterpolatorPtr_
     (
         master_
     ?
     new PrimitivePatchInterpolation
     <
             PrimitivePatch<face, Foam::List, pointField>
     >(*masterFaceZonePatchPtr_)
     :
     NULL
     ),
     slaveFaceZonePatchInterpolatorPtr_
     (
         master_
     ?
         new PrimitivePatchInterpolation
     <
             PrimitivePatch<face, Foam::List, pointField>
     >(*slaveFaceZonePatchPtr_)
     :
     NULL
     ),
    oldMasterFaceZonePoints_
    (
     master_ ?
     masterFaceZonePatchPtr_->localPoints()
     :
     vectorField(0,vector::zero)
     ),
    oldSlaveFaceZonePoints_
    (
     master_ ?
     slaveFaceZonePatchPtr_->localPoints()
     :
     vectorField(0,vector::zero)
     ),
    alg_
    (
        master_
        ? intersection::algorithmNames_.read
        (
            dict.lookup("projectionAlgo")
            ) : Foam::intersection::VISIBLE
        ),
    dir_(
        master_
        ? intersection::directionNames_.read
        (
            dict.lookup("projectionDir")
            ) : Foam::intersection::CONTACT_SPHERE
        ),
    curTimeIndex_(-1),
    iCorr_(0),
    correctionFreq_(master_ ? readInt(dict.lookup("correctionFrequency")) : 1),
    orthotropic_(false),
    stickSlipFieldPtr_
    (
     this->db().objectRegistry::foundObject<volScalarField>("stickSlipRegions")
     ?
     &this->db().objectRegistry::lookupObject<volScalarField>
     ("stickSlipRegions")
     :
     new volScalarField
     (
      IOobject
      (
       "stickSlipRegions",
       patch().boundaryMesh().mesh().time().timeName(),
       patch().boundaryMesh().mesh(),
       IOobject::READ_IF_PRESENT,
       IOobject::AUTO_WRITE
       ),
      patch().boundaryMesh().mesh(),
      dimensionedScalar("zero", dimless, 0)
      )
     ),
    forceCorrection_(false),
    nonLinear_(nonLinearGeometry::OFF)
{
    Info << "Creating " << solidContactFvPatchVectorField::typeName << " patch"
        << endl;

    // check shadow patch exists
    if (shadowPatchID_ == -1)
    {
        FatalError
            << "\nCannot find shadowPatch called " << dict.lookup("shadowPatch")
            << " for patch " << patch().name() << exit(FatalError);
    }

    if (dict.found("orthotropic"))
    {
        orthotropic_ = Switch(dict.lookup("orthotropic"));
        Info << "\torthotropic: " << orthotropic_ << endl;
    }

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

  // master creates contact laws
  // and the slave grabs the contact law pointers in updateCoeffs
    if (master_)
    {
        rigidMaster_ = Switch(dict.lookup("rigidMaster"));

        if (contactActive_)
        {
            const fvMesh& mesh = patch().boundaryMesh().mesh();
            label masterFaceZoneID =
                mesh.faceZones().findZoneID(masterFaceZoneName_);
            label slaveFaceZoneID =
                mesh.faceZones().findZoneID(slaveFaceZoneName_);
            if (masterFaceZoneID == -1)
            {
                FatalError
                    << "face zone " << masterFaceZoneName_ << " not found!"
                    << exit(FatalError);
            }
            else if (slaveFaceZoneID == -1)
            {
                FatalError
                    << "face zone " << masterFaceZoneName_ << " not found!"
                    << exit(FatalError);
            }

            normalContactModelPtr_ =
                normalContactModel::New
                (
                    dict.lookup("normalContactModel"),
                    patch(),
                    dict,
                    patch().index(), // master
                    shadowPatchID_, // slave
                    masterFaceZoneID_, // master face zone ID
                    slaveFaceZoneID_, // slave face zone ID
                    *masterFaceZonePatchPtr_,
                    *slaveFaceZonePatchPtr_
                    ).ptr();
            frictionContactModelPtr_ =
                frictionContactModel::New
                (
                    dict.lookup("frictionContactModel"),
                    patch(),
                    dict,
                    patch().index(), // master
                    shadowPatchID_, // slave
                    masterFaceZoneID_, // master face zone ID
                    slaveFaceZoneID_ // slave face zone ID
                    ).ptr();

            //- create zoneToZone or ggiZone for interpolation of traction from
            //  slave to master only required if the contact is not rigid
            if (interpolationMethod_ == "patchToPatch")
            {
                Info<< "\tInterpolation of traction from slave to master:"
                    << " patchToPatch"<< endl;
                // slaveToMasterPatchToPatchInterpolatorPtr_ =
                //   new zoneToZoneInterpolation
                //   (
                //    // slaveFaceZonePatch, // from zone
                //    // masterFaceZonePatch, // to zone
                //    *slaveFaceZonePatchPtr_, // from zone
                //    *masterFaceZonePatchPtr_, // to zone
                //    alg_,
                //    dir_
                //    );
                slaveToMasterPatchToPatchInterpolatorPtr_ =
                    new PatchToPatchInterpolation
                    <
                        PrimitivePatch<face, Foam::List, pointField>,
                        PrimitivePatch<face, Foam::List, pointField>
                        >
                    (
                        *slaveFaceZonePatchPtr_, // from zone
                        *masterFaceZonePatchPtr_, // to zone
                        alg_,
                        dir_
                    );
            }
            else if (interpolationMethod_ == "ggi")
            {
                // still created if master is rigid because it is used to
                // interpolate normals
                Info<< "\tInterpolation of traction from slave to master:"
                    << " ggi" << endl;
                slaveToMasterGgiInterpolatorPtr_ =
                    // new ggiZoneInterpolation
                    new GGIInterpolation<
                        PrimitivePatch<
                            face, Foam::List, pointField
                            >, PrimitivePatch< face, Foam::List, pointField > >
                    (
                        // masterFaceZonePatch, // master zone
                        // slaveFaceZonePatch, // slave zone
                        *masterFaceZonePatchPtr_, // master zone
                        *slaveFaceZonePatchPtr_, // slave zone
                        tensorField(0),
                        tensorField(0),
                        vectorField(0),
                        0.0,
                        0.0,
                        true,
                        ggiInterpolation::AABB
                    );
            }
            else
            {
                SeriousError
                    << "\n\nTraction interpolation method '"
                    << interpolationMethod_ << "' not found!\n"
                    << "\nValid interpolation methods are:\n"
                    << "ggi\n"
                    << "patchToPatch"
                    << endl
                    << exit(FatalError);
            }
        }
        else
        {
            Info<< "\t**The contact is not active, contact patches treated "
                << "as traction-free**" << endl;
            valueFraction() = symmTensor::zero;
        }
    } // end of if master

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
      //vectorField n = patch().nf();
      this->valueFraction() = symmTensor::zero; //sqr(n);
    }

    if (dict.found("value"))
    {
        Field<vector>::operator=(vectorField("value", dict, p.size()));
    }
    else
    {
        Field<vector> normalValue = transform(valueFraction(), refValue());

        Field<vector> gradValue =
            this->patchInternalField() + refGrad()/this->patch().deltaCoeffs();

        Field<vector> transformGradValue =
            transform(I - valueFraction(), gradValue);

        Field<vector>::operator=(normalValue + transformGradValue);
    }
}


solidContactFvPatchVectorField::solidContactFvPatchVectorField
(
    const solidContactFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    directionMixedFvPatchVectorField(ptf, iF),
    fieldName_(ptf.fieldName_),
    master_(ptf.master_),
    contactActive_(ptf.contactActive_),
    rigidMaster_(ptf.rigidMaster_),
    normalContactModelPtr_(ptf.normalContactModelPtr_),
    frictionContactModelPtr_(ptf.frictionContactModelPtr_),
    shadowPatchID_(ptf.shadowPatchID_),
    patchInChargeOfCorrection_(ptf.patchInChargeOfCorrection_),
    masterFaceZoneName_(ptf.masterFaceZoneName_),
    slaveFaceZoneName_(ptf.slaveFaceZoneName_),
    masterFaceZoneID_(ptf.masterFaceZoneID_),
    slaveFaceZoneID_(ptf.slaveFaceZoneID_),
    // masterFaceZonePatchPoints_(ptf.masterFaceZonePatchPoints_),
    // slaveFaceZonePatchPoints_(ptf.slaveFaceZonePatchPoints_),
    // masterFaceZonePatchFaces_(ptf.masterFaceZonePatchFaces_),
    // slaveFaceZonePatchFaces_(ptf.slaveFaceZonePatchFaces_),
    masterFaceZonePatchPtr_(ptf.masterFaceZonePatchPtr_),
    slaveFaceZonePatchPtr_(ptf.slaveFaceZonePatchPtr_),
    interpolationMethod_(ptf.interpolationMethod_),
    slaveToMasterPatchToPatchInterpolatorPtr_
    (ptf.slaveToMasterPatchToPatchInterpolatorPtr_),
    slaveToMasterGgiInterpolatorPtr_(ptf.slaveToMasterGgiInterpolatorPtr_),
    masterFaceZonePatchInterpolatorPtr_
    (ptf.masterFaceZonePatchInterpolatorPtr_),
    slaveFaceZonePatchInterpolatorPtr_(ptf.slaveFaceZonePatchInterpolatorPtr_),
    oldMasterFaceZonePoints_(ptf.oldMasterFaceZonePoints_),
    oldSlaveFaceZonePoints_(ptf.oldSlaveFaceZonePoints_),
    alg_(ptf.alg_),
    dir_(ptf.dir_),
    curTimeIndex_(ptf.curTimeIndex_),
    iCorr_(ptf.iCorr_),
    correctionFreq_(ptf.correctionFreq_),
    orthotropic_(ptf.orthotropic_),
    stickSlipFieldPtr_(ptf.stickSlipFieldPtr_),
    forceCorrection_(ptf.forceCorrection_),
    nonLinear_(ptf.nonLinear_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Map from self
void solidContactFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
//     if (normalContactModelPtr_ == NULL)
//     {
//         FatalErrorIn("solidContactFvPatchVectorField::autoMap")
//             << "NULL contact normal law"
//             << abort(FatalError);
//     }

//     if (frictionContactModelPtr_ == NULL)
//     {
//         FatalErrorIn("solidContactFvPatchVectorField::autoMap")
//             << "NULL contact friction law"
//             << abort(FatalError);
//     }

    directionMixedFvPatchVectorField::autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void solidContactFvPatchVectorField::rmap
(
    const fvPatchField<vector>& ptf,
    const labelList& addr
)
{
    directionMixedFvPatchVectorField::rmap(ptf, addr);

    // not sure if pointers are mapped correctly
    // be careful when there are topological changes to the patch
}


void solidContactFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (contactActive_)
    {
        // for the first updateCoeffs, the slave needs to grab the conatct
        // law pointers
        if (!normalContactModelPtr_)
        {
            if (!master_)
            {
                // grab pointers
                const volVectorField& Ufield =
                    this->db().objectRegistry::lookupObject<volVectorField>
                    (
                        fieldName_
                        );

                const solidContactFvPatchVectorField& Upatch =
                    refCast<const solidContactFvPatchVectorField>
                    (
                        Ufield.boundaryField()[shadowPatchID_]
                        );
                Info<< "\tSlave contact patch " << patch().name()
                    << " grabbing normalContactModel pointer from master"
                    << endl;
                normalContactModelPtr_ = Upatch.normalContactModelPtr();
                if (!normalContactModelPtr_)
                {
                    FatalError
                        << "\nThe patch " << patch().name()
                        << " has a NULL normalContactModel pointer!"
                        << exit(FatalError);
                }

                Info<< "\tSlave contact patch " << patch().name()
                    << " grabbing frictionContactModel pointer from master"
                    << endl;
                frictionContactModelPtr_ = Upatch.frictionContactModelPtr();
                if (!frictionContactModelPtr_)
                {
                    FatalError
                        << "\nThe patch " << patch().name()
                        << " has a NULL frictionContactModel pointer!"
                        << exit(FatalError);
                }

                // lookup master's correction frequency
                correctionFreq_ = Upatch.correctionFreq();

                // get face zone patch pointers
                masterFaceZonePatchPtr_ = Upatch.masterFaceZonePatchPtr();
                slaveFaceZonePatchPtr_ = Upatch.slaveFaceZonePatchPtr();

                // get patchToPatch and GGI interpolator pointers
                // one of them is NULL
                slaveToMasterGgiInterpolatorPtr_ =
                    Upatch.slaveToMasterGgiInterpolatorPtr();
                slaveToMasterPatchToPatchInterpolatorPtr_ =
                    Upatch.slaveToMasterPatchToPatchInterpolatorPtr();
            }
            else
            {
                FatalErrorIn("solidContactFvPatchVectorField::updateCoeff()")
                    << "NULL contactLaw\ncontactLaw not created by master."
                    << abort(FatalError);
            }
        }

        // if it is a new time step then reset iCorr
        if (curTimeIndex_ != this->db().time().timeIndex())
        {
            curTimeIndex_ = this->db().time().timeIndex();
            iCorr_ = 0;

            // update old face zone points
            if (master_ && nonLinear_ == nonLinearGeometry::UPDATED_LAGRANGIAN)
            {
                oldMasterFaceZonePoints_ =
                    masterFaceZonePatchPtr_->localPoints();
                oldSlaveFaceZonePoints_ =
                    slaveFaceZonePatchPtr_->localPoints();
            }
        }

        // the time taken to correct the contact may not be negligible
        // so reduce the correctiion frequency can speed up the simulation
        if (iCorr_++ % correctionFreq_ == 0
            || forceCorrection_ )
        {
            forceCorrection_ = false;

            // master moves faceZone patches to the current deformed position
            if (master_)
            {
                moveFaceZonePatches();
            }

            // just one of the patches updates the interpolators and
            // corrects the contact laws
            if (master_)
            {
                // update interpolation classes
                if (!rigidMaster_)
                {
                    if (slaveToMasterGgiInterpolatorPtr_)
                    {
                        slaveToMasterGgiInterpolatorPtr_->movePoints
                        (
                            tensorField(0),
                            tensorField(0),
                            vectorField(0)
                        );
                    }
                    else if (slaveToMasterPatchToPatchInterpolatorPtr_)
                    {
                        slaveToMasterPatchToPatchInterpolatorPtr_->movePoints();
                    }
                    else
                    {
                        FatalErrorIn("solidContactFvPatchVectorField::"
                                     "updateCoeff()")
                            << "Neither patchToPatch or GGI interpolator found!"
                            << abort(FatalError);
                    }
                }

                // correct laws
                // the normal model sets what face normal to use
                // on the slave e.g. use actual face normals, or master normals
                // interpolated to the slave, undeformed/deformed normals.
                // the friction model then uses these face normals
                vectorField slaveFaceNormals
                    (
                        patch().boundaryMesh().mesh().boundary()
                        [shadowPatchID_].size(),
                        vector::zero
                        );

                normalContactModelPtr_->correct
                    (
                        *masterFaceZonePatchPtr_,
                        *slaveFaceZonePatchPtr_,
                        alg_,
                        dir_,
                        fieldName_,
                        orthotropic_,
                        nonLinearGeometry::nonLinearNames_[nonLinear_],
                        slaveFaceNormals,
                        slaveToMasterGgiInterpolatorPtr_
                        );

                frictionContactModelPtr_->correct
                    (
                        normalContactModelPtr_->slavePressure(),
                        *masterFaceZonePatchPtr_,
                        *slaveFaceZonePatchPtr_,
                        alg_,
                        dir_,
                        interpolationMethod_,
                        fieldName_,
                        orthotropic_,
                        nonLinearGeometry::nonLinearNames_[nonLinear_],
                        slaveFaceNormals
                        );
            }

            // set boundary conditions
            bool incremental(fieldName_ == "DU");

            if (!master_)
            {
                // set refValue, refGrad and valueFraction

                // refValue
                refValue() =
                    normalContactModelPtr_->slaveDisp()
                    + frictionContactModelPtr_->slaveDisp();

                //refGrad - set traction
                refGrad() =
                    tractionBoundaryGradient::snGrad
                    (
                        frictionContactModelPtr_->slaveTraction()
                        + normalContactModelPtr_->slavePressure(),
                        scalarField(patch().size(),0.0),
                        fieldName_,
                        "U",
                        patch(),
                        orthotropic_,
                        nonLinear_,
                        incremental
                    );

                //valueFraction
                valueFraction() =
                    normalContactModelPtr_->slaveValueFrac()
                    + frictionContactModelPtr_->slaveValueFrac();
            }
            else
            {
                // master is always traction condition
                // interpolate traction from slave
                // Dirichlet-Neumann uses lagged traction
                // penalty uses same as for slave traction

                if (rigidMaster_)
                {
                    // set to master to traction free if it is rigid
                    refGrad() =
                        tractionBoundaryGradient::snGrad
                        (
                            vectorField(patch().size(),vector::zero),
                            scalarField(patch().size(),0.0),
                            fieldName_,
                            "U",
                            patch(),
                            orthotropic_,
                            nonLinear_,
                            incremental
                        );
                }
                else
                {
                    refGrad() = tractionBoundaryGradient::snGrad
                    (
                        interpolateSlaveToMaster
                        (
                            -frictionContactModelPtr_->slaveTraction()
                            -normalContactModelPtr_->slavePressure()
                        ),
                        scalarField(patch().size(),0.0),
                        fieldName_,
                        "U",
                        patch(),
                        orthotropic_,
                        nonLinear_,
                        incremental
                    );
                }
            }
        } // if correction freqeuncy
    } // if contactActive

    // if the contact is not active then the patches behave as traction free
    else
    {
        // set refGrad to traction free for master and slave
        bool incremental(fieldName_ == "DU");

        refGrad() =
            tractionBoundaryGradient::snGrad
            (
                vectorField(patch().size(),vector::zero),
                scalarField(patch().size(),0.0),
                fieldName_,
                "U",
                patch(),
                orthotropic_,
                nonLinear_,
                incremental
            );
    }

    directionMixedFvPatchVectorField::updateCoeffs();
}


// Interpolate traction from slave to master
tmp<vectorField> solidContactFvPatchVectorField::interpolateSlaveToMaster
(
 const vectorField slaveField
)
{
    if (!master_)
    {
      FatalError
          << "only the master may call the function "
          "solidContactFvPatchVectorField::interpolateSlaveToMaster"
          << exit(FatalError);
    }

  const fvMesh& mesh = patch().boundaryMesh().mesh();

  vectorField masterField
      (
          mesh.boundaryMesh()[patch().index()].size(),
          vector::zero
          );

  // for now, the slave must be the slave
  const label slaveStart
    = mesh.boundaryMesh()[shadowPatchID_].start();

  // global slave field
  vectorField globalSlaveField(slaveFaceZonePatchPtr_->size(), vector::zero);

  // put local slaveField into globalSlaveField
  forAll(slaveField, i)
    {
      globalSlaveField[
          mesh.faceZones()[slaveFaceZoneID_].whichFace(slaveStart + i)
          ] =
        slaveField[i];
    }

  //- exchange parallel data
  // sum because each face is only on one proc
  reduce(globalSlaveField, sumOp<vectorField>());

  // select interpolator - patchToPatch or GGI
  vectorField globalMasterInterpField
      (
          masterFaceZonePatchPtr_->size(),
          vector::zero
      );
  if (slaveToMasterPatchToPatchInterpolatorPtr_)
  {
      globalMasterInterpField =
          slaveToMasterPatchToPatchInterpolatorPtr_->faceInterpolate<vector>
          (
              globalSlaveField
          );
  }
  else if (slaveToMasterGgiInterpolatorPtr_)
  {
      globalMasterInterpField =
          slaveToMasterGgiInterpolatorPtr_->slaveToMaster
          (
              globalSlaveField
          );
  }
  else
  {
      FatalErrorIn("solidContactFvPatchVectorField::interpolateSlaveToMaster()")
          << "interpolationMethod is not patchToPatch or GGI!"
          << abort(FatalError);
  }

  // now put global back into local
  const label masterPatchStart
    = mesh.boundaryMesh()[patch().index()].start();

  tmp<vectorField> tmasterInterpField
      (
          new vectorField(masterFaceZonePatchPtr_->size(),vector::zero)
          );
  vectorField& masterInterpField = tmasterInterpField();

  forAll(masterInterpField, i)
    {
      masterInterpField[i] =
        globalMasterInterpField
        [
         mesh.faceZones()[masterFaceZoneID_].whichFace(masterPatchStart + i)
         ];
      }

    return tmasterInterpField;
}


//  Move the contact face zone patches to the deformed position
void solidContactFvPatchVectorField::moveFaceZonePatches()
{
      //OFstream outFile("localFaces_"+name(Pstream::myProcNo()));
      //outFile << slaveFaceZonePatchPtr_->localFaces() << endl;

  // method: we get the total displacement field for the global
  // face zone patches. We then interpolate these face values
  // to the vertices. And we move the vertices by these
  // interpolated displacements, so the global face zone patches
  // should be in the same deformed position on all procs.
  if (!master_)
    {
      FatalError
          << "Only the master may call the function "
          "solidContactFvPatchVectorField::moveFaceZonePatches"
          << exit(FatalError);
    }

  // update face zone patch interpolators
  masterFaceZonePatchInterpolatorPtr_->movePoints();
  slaveFaceZonePatchInterpolatorPtr_->movePoints();

  // get total displacement fields for the master and slave face zone patches
  const fvMesh& mesh = patch().boundaryMesh().mesh();
  vectorField globalMasterU(masterFaceZonePatchPtr_->size(), vector::zero);
  vectorField globalSlaveU(slaveFaceZonePatchPtr_->size(), vector::zero);
  const label masterPatchStart
    = mesh.boundaryMesh()[patch().index()].start();
  const label slavePatchStart
    = mesh.boundaryMesh()[shadowPatchID_].start();

  // get local total displacement fields
  const volVectorField& dispField =
    this->db().objectRegistry::lookupObject<volVectorField>(fieldName_);
  vectorField localMasterU = dispField.boundaryField()[patch().index()];
  vectorField localSlaveU = dispField.boundaryField()[shadowPatchID_];
  if (fieldName_ == "DU"
     &&
     nonLinear_ != nonLinearGeometry::UPDATED_LAGRANGIAN)
    {
      const volVectorField& totalDispField =
        this->db().objectRegistry::lookupObject<volVectorField>("U");
      localMasterU += totalDispField.boundaryField()[patch().index()];
      localSlaveU += totalDispField.boundaryField()[shadowPatchID_];
    }
  else if (fieldName_ != "U" && fieldName_ != "DU")
    {
      FatalError << "Displacement field must be U or DU!"
                 << exit(FatalError);
    }

  // add on initial position so that localSlaveU becomes current local
  // deformed position - done after
  //localMasterU += mesh.boundaryMesh()[patch().index()].faceCentres();
  //localSlaveU += mesh.boundaryMesh()[shadowPatchID_].faceCentres();

  // put localMasterU field into globalMasterU
  forAll(localMasterU, i)
    {
      globalMasterU[
          mesh.faceZones()[masterFaceZoneID_].whichFace(masterPatchStart + i)
          ] = localMasterU[i];
    }
  // put localSlaveU field into globalSlaveU
  forAll(localSlaveU, i)
    {
      globalSlaveU[
          mesh.faceZones()[slaveFaceZoneID_].whichFace(slavePatchStart + i)
          ] = localSlaveU[i];
    }

  //- exchange parallel data
  // sum because each face is only on one proc
  reduce(globalMasterU, sumOp<vectorField>());
  reduce(globalSlaveU, sumOp<vectorField>());

  // hmnnn - globalSlaveU is exactly the same on every proc
  // but slaveFaceZonePatchInterpolatorPtr_ needs globalSlaveU to
  // be ordered in locally to be in the same order as the local
  // globalFaceZone

  // interpolate displacement from face centre to vertices
  vectorField globalMasterNewPoints =
    masterFaceZonePatchInterpolatorPtr_->faceToPointInterpolate(globalMasterU);
  vectorField globalSlaveNewPoints =
    slaveFaceZonePatchInterpolatorPtr_->faceToPointInterpolate(globalSlaveU);

  // debug: try x motion to zero
  // vector iHat(1,0,0);
  // globalSlaveNewPoints = (I - sqr(iHat)) & globalSlaveNewPoints;
  // globalMasterNewPoints = (I - sqr(iHat)) & globalMasterNewPoints;

  // Add displacement to undeformed mesh points
  // globalMasterNewPoints +=
  // mesh.faceZones()[masterFaceZoneID_]().localPoints();
  // globalSlaveNewPoints += mesh.faceZones()[slaveFaceZoneID_]().localPoints();
  // BUG FIX 13-08-13 - philipc
  // the faceZones keep with the mesh (mesh.faceZones()) are not moved correctly
  // so I must keep a copy of the old time-step points
  globalMasterNewPoints += oldMasterFaceZonePoints_;
  globalSlaveNewPoints += oldSlaveFaceZonePoints_;

  // find points which are on symmetryPlanes and force them
  // to stay exactly on the symmetryPlane - WIP
  // const labelList& slaveBoundaryPoints =
  // slaveFaceZonePatchPtr_->boundaryPoints();
  // const labelList& slaveMeshPoints = slaveFaceZonePatchPtr_->meshPoints();
  // labelList slaveBoundaryPointsGlobalIndex(slaveBoundaryPoints.size(), -1);
  // forAll(slaveBoundaryPointsGlobalIndex, pointi)
  //   {
  //     slaveBoundaryPointsGlobalIndex[pointi] =
  // slaveMeshPoints[slaveBoundaryPoints[pointi]];
  //   }

  // as I can't figure out how to move the points
  // of the current face zone patches
  // I will delete the face zone patches and create new ones

  // delete old face zone patches
  delete masterFaceZonePatchPtr_;
  delete slaveFaceZonePatchPtr_;

  // create new face zone patches with deformed points
  masterFaceZonePatchPtr_ =
    new PrimitivePatch<face, Foam::List, pointField>
    (
     mesh.faceZones()[masterFaceZoneID_]().localFaces(),
     globalMasterNewPoints
     );

  // masterFaceZonePatchPtr_->writeVTK
  //   (
  //         fileName("masterFaceZonePatch"+name(Pstream::myProcNo())),
  //    masterFaceZonePatchPtr_->localFaces(),
  //    masterFaceZonePatchPtr_->localPoints()
  //    );

  slaveFaceZonePatchPtr_ =
    new PrimitivePatch<face, Foam::List, pointField>
    (
     mesh.faceZones()[slaveFaceZoneID_]().localFaces(),
     globalSlaveNewPoints
     );

  // slaveFaceZonePatchPtr_->writeVTK
  //   (
  //    fileName("slaveFaceZonePatch"+name(Pstream::myProcNo())),
  //    slaveFaceZonePatchPtr_->localFaces(),
  //    slaveFaceZonePatchPtr_->localPoints()
  //    );

  // The patchToPatch, GGI and primitivePatch interpolators point to the
  // old primitive patches so I must delete these interpolators
  // it would be much nicer if I could just move the faceZonePatch points...
  if (slaveToMasterPatchToPatchInterpolatorPtr_)
  {
      delete slaveToMasterPatchToPatchInterpolatorPtr_;
      slaveToMasterPatchToPatchInterpolatorPtr_ =
          new PatchToPatchInterpolation<
              PrimitivePatch<
                  face, Foam::List, pointField
                  >, PrimitivePatch<face, Foam::List, pointField> >
          (
              *slaveFaceZonePatchPtr_, // from zone
              *masterFaceZonePatchPtr_, // to zone
              alg_,
              dir_
          );
  }
  else if (slaveToMasterGgiInterpolatorPtr_)
  {
      delete slaveToMasterGgiInterpolatorPtr_;
      slaveToMasterGgiInterpolatorPtr_ =
        new GGIInterpolation<
            PrimitivePatch<
                face, Foam::List, pointField
                >, PrimitivePatch< face, Foam::List, pointField > >
          (
              *masterFaceZonePatchPtr_, // master zone
              *slaveFaceZonePatchPtr_, // slave zone
              tensorField(0),
              tensorField(0),
              vectorField(0),
              0.0,
              0.0,
              true,
              ggiInterpolation::AABB
          );
  }

  // and primitive patch interpolators
  delete masterFaceZonePatchInterpolatorPtr_;
  masterFaceZonePatchInterpolatorPtr_ =
    new PrimitivePatchInterpolation<
        PrimitivePatch<face, Foam::List, pointField>
        >(*masterFaceZonePatchPtr_);
  delete slaveFaceZonePatchInterpolatorPtr_;
  slaveFaceZonePatchInterpolatorPtr_ =
    new PrimitivePatchInterpolation<
        PrimitivePatch<face, Foam::List, pointField>
        >(*slaveFaceZonePatchPtr_);

  // Also maybe I should correct motion for 2D models
  // OK for now
}


// check shadow patch and face zones exist
bool solidContactFvPatchVectorField::checkPatchAndFaceZones
(const dictionary& dict) const
{
  // check shadow patch
  word shadowName = dict.lookup("shadowPatch");
  label shadowPatchID = patch().patch().boundaryMesh().findPatchID(shadowName);
  if (shadowPatchID == -1)
    {
      FatalError << "shadowPatch " << shadowName << " not found for patch "
                 << patch().name() << exit(FatalError);
    }

  word curZoneName = patch().name()+"FaceZone";
  word shadowZoneName =
  patch().boundaryMesh().mesh().boundary()[shadowPatchID].name() + "FaceZone";
  label curPatchFaceZoneID =
    patch().boundaryMesh().mesh().faceZones().findZoneID(curZoneName);
  if (curPatchFaceZoneID == -1)
    {
      FatalError
          << "faceZone " << curZoneName
          << " not found and is required for the solidContact boundaries\n"
          << "To create a faceZone from a patch, use the setSet and "
          << "setsToZone utilities:\n"
          << "\tsetSet\n"
          << "\tfaceSet "<<curZoneName<<" new patchToFace "<<patch().name()<< nl
          << "\tfaceSet "<<shadowZoneName<<" new patchToFace "<<shadowName<<nl
          << "\tquit\n"
          << "\tsetsToZones -noFlipMap\n"
          << exit(FatalError);
    }

  label shadowPatchFaceZoneID =
  patch().boundaryMesh().mesh().faceZones().findZoneID(shadowZoneName);

  if (shadowPatchFaceZoneID == -1)
  {
      FatalError
          << "faceZone " << shadowZoneName
          << " not found and is required for the solidContact boundaries\n"
          << "To create a faceZone from a patch, use the setSet and "
          << "setsToZone utilities:\n"
          << "\tsetSet\n"
          << "\tfaceSet "<<curZoneName<<" new patchToFace "<<patch().name()<< nl
          << "\tfaceSet "<<shadowZoneName<<" new patchToFace "<<shadowName<<nl
          << "\tquit\n"
          << "\tcreate zones from sets\n"
          << "\tsetsToZones -noFlipMap\n"
          << exit(FatalError);
  }

  // if the total amount of faces in the master or slave face zones is zero
  // then something is
  // wrong with the face zones - they were probably created wrong.
  if (
      returnReduce
      (
          patch().boundaryMesh().mesh().faceZones()[curPatchFaceZoneID].size(),
          sumOp<label>()
          ) < 1
      )
  {
      FatalError
          << "faceZone " << curZoneName
          << ", which is required for the solidContact boundaries,"
          << " has no faces!\n"
          << "You probably made a mistake creating the faceZones."
          << exit(FatalError);
  }

  if (
      returnReduce
      (
          patch().boundaryMesh().mesh().faceZones()
          [shadowPatchFaceZoneID].size(),
          sumOp<label>()
          ) < 1)
  {
      FatalError
          << "faceZone " << shadowZoneName
          << ", which is required for the solidContact boundaries,"
          << " has no faces!\n"
          << "You probably made a mistake creating the faceZones."
          << exit(FatalError);
  }

  return true;
}

void solidContactFvPatchVectorField::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    vectorField normalValue = transform(valueFraction(), refValue());

    //- non-orthogonal correction vectors needed to calculate gradValue
    const fvPatchField<tensor>& gradField =
        patch().lookupPatchField<volTensorField, tensor>
        (
            "grad(" + fieldName_ + ")"
        );
    vectorField n = patch().nf();
    vectorField delta = patch().delta();
    vectorField k = delta - n*(n&delta);

    vectorField gradValue =
      this->patchInternalField()
      + (k&gradField.patchInternalField())
      + refGrad()/this->patch().deltaCoeffs();

    Field<vector> transformGradValue =
      transform(I - valueFraction(), gradValue);

    Field<vector>::operator=(normalValue + transformGradValue);

    fvPatchField<vector>::evaluate();
}

Foam::tmp<Foam::Field<vector> > solidContactFvPatchVectorField::
snGrad() const
{
    vectorField pif = this->patchInternalField();

    vectorField normalValue = transform(valueFraction(), refValue());

    //- non-orthogonal correction vectors needed to calculate gradValue
    const fvPatchField<tensor>& gradField =
        patch().lookupPatchField<volTensorField, tensor>
        (
            "grad(" + fieldName_ + ")"
        );
    vectorField n = patch().nf();
    vectorField delta = patch().delta();
    vectorField k = delta - n*(n&delta);

    vectorField gradValue =
      this->patchInternalField()
      + (k&gradField.patchInternalField())
      + refGrad()/this->patch().deltaCoeffs();

    vectorField transformGradValue =
      transform(I - valueFraction(), gradValue);

    return
      (
       (normalValue + transformGradValue)
       - (pif + (k&gradField.patchInternalField()))
       )*this->patch().deltaCoeffs();
}


//- Increment of dissipated energy due to friction
tmp<scalarField> solidContactFvPatchVectorField::Qc() const
{
    tmp<scalarField> tQc(new scalarField(patch().size(), 0.0));
    scalarField& Qc = tQc();

    // Integrate energy using trapezoidal rule
    // 0.5*averageForce*incrementOfDisplacement

    // displacement increment
    const vectorField& curDU = *this;

    // average force
    const fvPatchField<tensor>& gradField =
        patch().lookupPatchField<volTensorField, tensor>
        (
            "grad(" + fieldName_ + ")"
        );

    // current Cauchy traction for nonlinear
    bool incremental(fieldName_ == "DU");

    const vectorField curTrac =
        tractionBoundaryGradient::traction
        (
            gradField,
            fieldName_,
            "U",
            patch(),
            orthotropic_,
            nonLinear_,
            incremental
        );

    vectorField Sf = patch().Sf();
    if
    (
        nonLinear_ != nonLinearGeometry::OFF
     && nonLinear_ == nonLinearGeometry::TOTAL_LAGRANGIAN
    )
    {
        // current areas
        const fvPatchField<tensor>& gradU =
            patch().lookupPatchField<volTensorField, tensor>("grad(U)");
        tensorField F = I + gradField + gradU;
        tensorField Finv = inv(F);
        scalarField J = det(F);
        Sf = J*(Finv & Sf);
    }
    const scalarField magSf = mag(Sf);

    const vectorField curForce = magSf*curTrac;

    const fvPatchField<symmTensor>& oldSigma =
        patch().lookupPatchField<volSymmTensorField, symmTensor>("sigma_0");

    vectorField oldSf = patch().Sf();
    if (nonLinear_ != nonLinearGeometry::OFF)
    {
        // current areas
        //tensorField F = I + gradField;
        if (nonLinear_ == nonLinearGeometry::TOTAL_LAGRANGIAN)
        {
            const fvPatchField<tensor>& gradU =
                patch().lookupPatchField<volTensorField, tensor>("grad(U)");
            tensorField F = I + gradU;
            tensorField Finv = inv(F);
            scalarField J = det(F);
            // rotate initial reference area to area of previous time-step
            oldSf = J*(Finv & Sf);
        }
        else if (nonLinear_ == nonLinearGeometry::UPDATED_LAGRANGIAN)
        {
            tensorField F = I + gradField;
            tensorField Finv = inv(F);
            scalarField J = det(F);
            // rotate current area to area of previous time-step
            oldSf = (F & Sf)/J;
        }
        else
        {
            FatalError << "solidContact::Qc() nonLinear type not known!"
                << exit(FatalError);
        }
    }

    const vectorField oldForce = oldSf & oldSigma;
    const vectorField avForce = 0.5*(oldForce + curForce);

    // Increment of dissiptaed frictional energy for this timestep
    Qc = mag(0.5*(avForce & curDU));

    return tQc;
}


// Write
void solidContactFvPatchVectorField::write(Ostream& os) const
{
  //Info << "writing..."<<flush;
    directionMixedFvPatchVectorField::write(os);

    os.writeKeyword("master") << master_ << token::END_STATEMENT << nl;
    os.writeKeyword("shadowPatch")
        << patch().boundaryMesh().mesh().boundary()[shadowPatchID_].name()
        << token::END_STATEMENT << nl;
    os.writeKeyword("orthotropic") << orthotropic_
                                   << token::END_STATEMENT << nl;
    os.writeKeyword("nonLinear")
        << nonLinearGeometry::nonLinearNames_[nonLinear_]
        << token::END_STATEMENT << nl;
    os.writeKeyword("contactActive")
        << contactActive_ << token::END_STATEMENT << nl;

    if (master_)
      {
        os.writeKeyword("rigidMaster")
            << rigidMaster_ << token::END_STATEMENT << nl;
        if (contactActive_)
          {
            os.writeKeyword("normalContactModel")
                << normalContactModelPtr_->type()
                << token::END_STATEMENT << nl;
            if (normalContactModelPtr_)
            {
                normalContactModelPtr_->writeDict(os);
            }
            os.writeKeyword("frictionContactModel")
                << frictionContactModelPtr_->type()
                << token::END_STATEMENT << nl;
            if (frictionContactModelPtr_)
            {
                frictionContactModelPtr_->writeDict(os);
            }
            os.writeKeyword("interpolationMethod")
                << interpolationMethod_ << token::END_STATEMENT << nl;
            os.writeKeyword("projectionAlgo")
                << intersection::algorithmNames_[alg_]
                << token::END_STATEMENT << nl;
            os.writeKeyword("projectionDir")
                << intersection::directionNames_[dir_]
                << token::END_STATEMENT << nl;
            os.writeKeyword("correctionFrequency")
                << correctionFreq_ << token::END_STATEMENT << nl;
          }
      }

    // set stickSlipFieldPtr_, it should be automatically written
    if (!master_ && contactActive_)
      {
        // apologies - I couldn't think of a tidy way to
        // pass around this field so I cast away the const-ness
        volScalarField* unConstStickSlipFieldPtr_ =
          const_cast<volScalarField*>(stickSlipFieldPtr_);

        if (!frictionContactModelPtr_)
          {
            // pointer gets dropped sometimes, not sure what the problem is
            // FatalError << "fricContactModel is NULL for slave"
            //              << exit(FatalError);
            Warning << "fricContactModel is NULL for slave" << nl
                    << "stickSlip field not written"
                    << endl;
          }
        else
          {
            unConstStickSlipFieldPtr_->boundaryField()[patch().index()] =
              frictionContactModelPtr_->stickSlipFaces();
          }
      }

    //Info << "done" << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, solidContactFvPatchVectorField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
