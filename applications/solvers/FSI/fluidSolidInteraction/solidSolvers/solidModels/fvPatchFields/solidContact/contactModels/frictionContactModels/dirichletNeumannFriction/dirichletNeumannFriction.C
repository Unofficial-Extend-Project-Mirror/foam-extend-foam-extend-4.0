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

#include "dirichletNeumannFriction.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include "PatchToPatchInterpolation.H"
#include "ggiInterpolation.H"
#include "tractionBoundaryGradient.H"
#include "fvc.H"
#include "primitivePatchInterpolation.H"
#include "nonLinearGeometry.H"

#include "solidContactFvPatchVectorField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dirichletNeumannFriction, 0);
    addToRunTimeSelectionTable
    (
        frictionContactModel,
        dirichletNeumannFriction,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::dirichletNeumannFriction::dirichletNeumannFriction
(
    const word& name,
    const fvPatch& patch,
    const dictionary& dict,
    const label masterPatchID,
    const label slavePatchID,
    const label masterFaceZoneID,
    const label slaveFaceZoneID
)
:
    frictionContactModel
    (
        name,
        patch,
        dict,
        masterPatchID,
        slavePatchID,
        masterFaceZoneID,
        slaveFaceZoneID
    ),
    frictionContactModelDict_(dict.subDict(name+"FrictionModelDict")),
    frictionLawPtr_(NULL),
    mesh_(patch.boundaryMesh().mesh()),
    slaveDisp_(mesh().boundaryMesh()[slavePatchID].size(), vector::zero),
    oldSlaveDisp_(slaveDisp_),
    slip_(slaveDisp_.size(), vector::zero),
    oldSlip_(slaveDisp_.size(), vector::zero),
    slipTol_
    (
        frictionContactModelDict_.lookupOrDefault<scalar>("slipTol", 1e-11)
    ),
    slaveTraction_(mesh().boundaryMesh()[slavePatchID].size(), vector::zero),
    oldSlaveTraction_(slaveTraction_),
    slaveTractionForMaster_(slaveTraction_.size(), vector::zero),
    oldSlaveTractionForMaster_(slaveTraction_.size(), vector::zero),
    slaveValueFrac_
    (
        mesh_.boundaryMesh()[slavePatchID].size(), symmTensor::zero
    ),
    oldSlaveValueFrac_(slaveValueFrac_),
    urDisp_
    (
        frictionContactModelDict_.lookupOrDefault<scalar>
        (
            "relaxationFactor", 0.1
        )
    ),
    urTrac_
    (
        frictionContactModelDict_.lookupOrDefault<scalar>
        (
            "relaxationFactorTraction", urDisp_
        )
    ),
    contactIterNum_(0),
    infoFreq_(readInt(frictionContactModelDict_.lookup("infoFrequency"))),
    oscillationCorr_(frictionContactModelDict_.lookup("oscillationCorrection")),
    smoothingSteps_
    (
        readInt(frictionContactModelDict_.lookup("smoothingSteps"))
    ),
    oldStickSlip_(slaveDisp_.size(), 0.0),
    contactFilePtr_(NULL)
{
    // Create friction law
    frictionLawPtr_ =
        frictionLaw::New
        (
            frictionContactModelDict_.lookup("frictionLaw"),
            *this,
            frictionContactModelDict_
        ).ptr();

    // Open contact info file
    if (Pstream::master())
    {
        word masterName = mesh_.boundary()[masterPatchID].name();
        word slaveName = mesh_.boundary()[slavePatchID].name();
        contactFilePtr_ =
            new OFstream
            (
                fileName("frictionContact_"+masterName+"_"+slaveName+".txt")
            );
        OFstream& contactFile = *contactFilePtr_;
        int width = 20;
        contactFile
            << "time";
        contactFile.width(width);
        contactFile
            << "iterNum";
        contactFile.width(width);
        contactFile
            << "relaxationFactor";
        contactFile.width(width);
        contactFile
            << "slipFaces";
        contactFile.width(width);
        contactFile
            << "stickFaces";
        contactFile.width(width);
        contactFile
            << "maxMagSlaveTraction" << endl;
    }
}


// Construct as a copy
Foam::dirichletNeumannFriction::dirichletNeumannFriction
(
    const dirichletNeumannFriction& fm
)
:
    frictionContactModel(fm),
    frictionContactModelDict_(fm.frictionContactModelDict_),
    frictionLawPtr_(fm.frictionLawPtr_->clone().ptr()),
    mesh_(fm.mesh_),
    slaveDisp_(fm.slaveDisp_),
    oldSlaveDisp_(fm.oldSlaveDisp_),
    slip_(fm.slip_),
    oldSlip_(fm.oldSlip_),
    slipTol_(fm.slipTol_),
    slaveTraction_(fm.slaveTraction_),
    oldSlaveTraction_(fm.oldSlaveTraction_),
    slaveTractionForMaster_(fm.slaveTractionForMaster_),
    oldSlaveTractionForMaster_(fm.oldSlaveTractionForMaster_),
    slaveValueFrac_(fm.slaveValueFrac_),
    oldSlaveValueFrac_(slaveValueFrac_),
    urDisp_(fm.urDisp_),
    urTrac_(fm.urTrac_),
    contactIterNum_(fm.contactIterNum_),
    infoFreq_(fm.infoFreq_),
    oscillationCorr_(fm.oscillationCorr_),
    smoothingSteps_(fm.smoothingSteps_),
    oldStickSlip_(fm.oldStickSlip_),
    contactFilePtr_(NULL)
{
    if (fm.contactFilePtr_)
    {
        contactFilePtr_ = new OFstream(*fm.contactFilePtr_);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::dirichletNeumannFriction::correct
(
    const vectorField& slavePressure,
    const primitiveFacePatch& masterFaceZonePatch,
    const primitiveFacePatch& slaveFaceZonePatch,
    const intersection::algorithm alg,
    const intersection::direction dir,
    const word interpolationMethod,
    const word fieldName,
    const Switch orthotropic,
    const nonLinearGeometry::nonLinearType nonLinear,
    const vectorField& slaveFaceNormals
)
{
    const fvMesh& mesh = mesh_;
    const label slavePatchIndex = slavePatchID();
    const label masterPatchIndex = masterPatchID();
    contactIterNum_++;

    // we have local masterDU and we want to interpolate it to the slave
    // to get local masterDUInterpToSlave (i.e. masterDU interpolated to the
    // slave)
    // so the method is:
    // create global masterDU field
    // interpolate global masterDU from master global face zone to slave global
    // zone then find local masterDUInterpToSlave from the global interpolated
    // field

    vectorField masterDUInterpToSlave
        (
            mesh.boundaryMesh()[slavePatchIndex].size(), vector::zero
        );

    // Global master DU
    vectorField globalMasterDU(masterFaceZonePatch.size(), vector::zero);

    // lookup current displacement field
    const volVectorField& dispField =
        mesh.objectRegistry::lookupObject<volVectorField>(fieldName);

    // local master and slave DU increment
    vectorField masterDU = dispField.boundaryField()[masterPatchIndex];
    vectorField slaveDU = dispField.boundaryField()[slavePatchIndex];

    if (fieldName == "U")
    {
        // lookup old U
        const volVectorField& dispOldField =
            mesh.objectRegistry::lookupObject<volVectorField>(fieldName + "_0");

        // subtract old U
        masterDU -= dispOldField.boundaryField()[masterPatchIndex];
        slaveDU -= dispOldField.boundaryField()[slavePatchIndex];
    }
    else if (fieldName != "DU")
    {
        FatalError
            << "dirichletNeumannFriction::correct()\n"
            " The displacement field must be called U or DU"
            << abort(FatalError);
    }

    // put local masterDU into globalMasterDU
    const label masterPatchStart
        = mesh.boundaryMesh()[masterPatchIndex].start();
    forAll(masterDU, i)
    {
        globalMasterDU
            [
                mesh.faceZones()[masterFaceZoneID()].whichFace
                    (
                        masterPatchStart + i
                    )
            ] = masterDU[i];
    }

    // exchange parallel data
    // sum because each face is only on one proc
    reduce(globalMasterDU, sumOp<vectorField>());

    // globalMasterDU is interpolated to the slave
    vectorField globalMasterDUInterpToSlave
        (
            slaveFaceZonePatch.size(), vector::zero
        );

    // interpolate DU from master to slave using inverse distance or ggi
    if (interpolationMethod == "patchToPatch")
    {
        PatchToPatchInterpolation<primitiveFacePatch, primitiveFacePatch>
            masterToSlavePatchToPatchInterpolator
            (
                masterFaceZonePatch, // from zone
                slaveFaceZonePatch, // to zone
                alg,
                dir
            );
        globalMasterDUInterpToSlave =
            masterToSlavePatchToPatchInterpolator.faceInterpolate<vector>
            (
                globalMasterDU
            );
    }
    else if (interpolationMethod == "ggi")
    {
        GGIInterpolation<primitiveFacePatch, primitiveFacePatch>
            masterToSlaveGgiInterpolator
            (
                masterFaceZonePatch, // master zone
                slaveFaceZonePatch, // slave zone
                tensorField(0),
                tensorField(0),
                vectorField(0),
                0.0,
                0.0,
                true,
                ggiInterpolation::AABB
            );
        globalMasterDUInterpToSlave =
            masterToSlaveGgiInterpolator.masterToSlave
            (
                globalMasterDU
            );
    }
    else
    {
        FatalError
            << "dirichletNeumannFriction::correct()" << nl
            << "interpolationMethod " << interpolationMethod << " not known"
            << nl << "interpolationMethod must be patchToPatch or ggi"
            << abort(FatalError);
    }

    // now put global back into local
    const label slavePatchStart
        = mesh.boundaryMesh()[slavePatchIndex].start();

    forAll(masterDUInterpToSlave, i)
    {
        masterDUInterpToSlave[i] =
            globalMasterDUInterpToSlave
            [
                mesh.faceZones()[slaveFaceZoneID()].whichFace
                (
                    slavePatchStart + i
                )
            ];
    }

    // Now masterDUInterpToSlave should have masterDU interpolated to the slave

    // Calculate current slave shear traction from the normal gradient field
    const fvPatch& slavePatch = mesh.boundary()[slavePatchIndex];
    const fvPatchField<tensor>& gradField =
        slavePatch.lookupPatchField<volTensorField, tensor>
            (
                "grad(" + fieldName + ")"
            );

    slaveTractionForMaster_ =
        (I - sqr(slaveFaceNormals))
        & tractionBoundaryGradient().traction
        (
            gradField,                 // grad field
            fieldName,                 // working field name
            "U",                       // total field name
            slavePatch,                // polyPatch
            bool(fieldName == "DU")    // incremental
        );


    // algorithm
    // if the face pressure is negative/zero then the friction is zero
    // so set valueFrac to zero and traction to zero
    // if the face pressure is positive and the shear traction is less
    // than fricCoeff*pressure then this is a sticking face so
    // set the valueFrac to (I-n^2) and set the disp to remove any
    // slip
    // if the face pressure is positive and the shear traction is greater
    // than fricCoeff*pressure then this is a slipping face so
    // set the valueFrac to zero and set the shear traction to
    // fricCoeff*pressure in the opposite direction to slip
    // if the shear traction on a slipping face is acting in the same
    // direction as the slip then this face should not be slipping
    // so we make it a sticking face
    // const volVectorField& prevSlaveDispField =
    //   mesh.objectRegistry::lookupObject<volVectorField>(fieldName);
    // const vectorField prevSlaveShearDisp =
    //     (I - sqr(slaveFaceNormals))
    //     & dispField.boundaryField()[slavePatchIndex];
    label numSlipFaces = 0;
    label numStickFaces = 0;
    scalarField& stickSlip = stickSlipFaces();
    //const scalar maxMagSlavePressure = gMax(magSlavePressure);
    // const scalarField magSlavePressure = -slaveFaceNormals & slavePressure;
    const scalarField magSlavePressure = mag(slavePressure);
    scalar maxMagSlavePressure = 0.0;
    if (slavePressure.size() > 0)
    {
        maxMagSlavePressure = max(magSlavePressure);
    }
    reduce(maxMagSlavePressure, maxOp<scalar>());

    // slip is the difference between the master tangential DU and slave
    // tangential DU
    slip_ = (I - sqr(slaveFaceNormals)) & (slaveDU - masterDUInterpToSlave);

    // Under-relax the slip
    //slip_ = urDisp_*slip_ + (1.0 - urDisp_)*oldSlip_;
    //oldSlip_ = slip_;

    // We only consider slip with a magnitude greater than the slipTol
    vectorField slipDir = slip_;
    scalarField magSlip = mag(slip_);
    forAll(slipDir, facei)
    {
        if (magSlip[facei] > slipTol_)
        {
            slipDir[facei] /= magSlip[facei];
        }
        else
        {
            slipDir[facei] = vector::zero;
            magSlip[facei] = 0.0;
        }
    }

    // const scalar pressureTol = max(1e-3*maxMagSlavePressure, 1e3);
    //const scalar pressureTol = max(1e-2*maxMagSlavePressure, 1e3);

    scalarField slipTrac = frictionLawPtr_->slipTraction(magSlavePressure);

    // lookup areaInContact
    const solidContactFvPatchVectorField& DUpatch =
        refCast<const solidContactFvPatchVectorField>
        (
            dispField.boundaryField()[slavePatchIndex]
        );
    const scalarField& areaInContact = DUpatch.normalModel().areaInContact();

    forAll(slaveDisp_, facei)
    {
        // Shear traction is zero outside of contact region; here we use the
        // slavePressure to decide where the contact area is. We may need to use
        // a better way
        // Try using areaInContact
        //if (mag(slavePressure[facei]) < pressureTol)
        if (areaInContact[facei] < SMALL)
        {
            // not in contact
            //slaveDisp_[facei] = prevSlaveShearDisp[facei];
            //slaveDisp_[facei] = vector::zero;
            slaveTraction_[facei] = vector::zero;
            slaveTractionForMaster_[facei] = vector::zero;
            slaveValueFrac_[facei] = symmTensor::zero;

            stickSlip[facei] = -1;
        }
        else if
        (
            (mag(slaveTractionForMaster_[facei]) > 0.999*slipTrac[facei])
            && // opposite directions
            ((slip_[facei] & slaveTractionForMaster_[facei]) < 0.0)
        )
        {
            // Slip
            //slaveDisp_[facei] = -slip_[facei] + prevSlaveShearDisp[facei];
            //slaveDisp_[facei] -= slip_[facei];
            slaveTraction_[facei] = -slipDir[facei]*slipTrac[facei];
            slaveValueFrac_[facei] = symmTensor::zero;
            numSlipFaces++;

            stickSlip[facei] = 0;
        }
        else
        {
            // stick
            //slaveDisp_[facei] = -slip_[facei] + prevSlaveShearDisp[facei];
            slaveDisp_[facei] -= slip_[facei];
            slaveTraction_[facei] = -slipDir[facei]*slipTrac[facei];
            //slaveTraction_[facei] = slaveShearTraction[facei];
            //slaveTraction_[facei] = vector::zero;
            slaveValueFrac_[facei] = I - sqr(slaveFaceNormals[facei]);
            numStickFaces++;

            stickSlip[facei] = 2;
        }
    }

    // Correct oscillations
    if (oscillationCorr_)
    {
        // interpolate face values to points then interpolate back
        // this essentially smooths the field
        primitivePatchInterpolation localSlaveInterpolator
            (
                mesh.boundaryMesh()[slavePatchIndex]
            );
        vectorField slavePointValues
            (
                mesh.boundaryMesh()[slavePatchIndex].nPoints(), vector::zero
            );

        for (int i = 0; i < smoothingSteps_; i++)
        {
            slavePointValues =
                localSlaveInterpolator.faceToPointInterpolate<vector>
                    (
                        slaveDisp_
                    );
            slaveDisp_ =
                localSlaveInterpolator.pointToFaceInterpolate<vector>
                    (
                        slavePointValues
                    );
            slavePointValues =
                localSlaveInterpolator.faceToPointInterpolate<vector>
                    (
                        slaveTraction_
                    );
            slaveTraction_ =
                localSlaveInterpolator.pointToFaceInterpolate<vector>
                    (
                        slavePointValues
                    );
            // slavePointValues =
            //     localSlaveInterpolator.faceToPointInterpolate<vector>
            //         (
            //             slaveTractionForMaster_
            //         );
            // slaveTractionForMaster_ =
            //     localSlaveInterpolator.pointToFaceInterpolate<vector>
            //         (
            //             slavePointValues
            //         );

            // make sure no normal component
            slaveDisp_ = (I - sqr(slaveFaceNormals)) & slaveDisp_;
            slaveTraction_ = (I - sqr(slaveFaceNormals)) & slaveTraction_;
            // slaveTractionForMaster_ =
            //     (I - sqr(slaveFaceNormals)) & slaveTractionForMaster_;
        }
    }

    // Under-relaxation
    //slaveDisp_ = urDisp_*slaveDisp_ + (1.0 - urDisp_)*prevSlaveShearDisp;
    slaveDisp_ = urDisp_*slaveDisp_ + (1.0 - urDisp_)*oldSlaveDisp_;
    oldSlaveDisp_ = slaveDisp_;
    slaveValueFrac_ =
        urDisp_*slaveValueFrac_ + (1.0 - urDisp_)*oldSlaveValueFrac_;
    oldSlaveValueFrac_ = slaveValueFrac_;
    slaveTraction_ = urTrac_*slaveTraction_ + (1.0 - urTrac_)*oldSlaveTraction_;
    oldSlaveTraction_ = slaveTraction_;
    slaveTractionForMaster_ =
        urTrac_*slaveTractionForMaster_
        + (1.0 - urTrac_)*oldSlaveTractionForMaster_;
    oldSlaveTractionForMaster_ = slaveTractionForMaster_;
    stickSlip = urDisp_*stickSlip + (1.0 - urDisp_)*oldStickSlip_;
    oldStickSlip_ = stickSlip;

    // get global values
    // in parallel, the log is poluted with warnings that
    // I am getting max of a list of size zero so
    // I will get the max of procs which have some
    // of the slave faces
    //scalar maxMagMasterTraction = gMax(mag(slaveTraction_))
    scalar maxMagMasterTraction = 0.0;
    if (slaveTraction_.size() > 0)
    {
        maxMagMasterTraction = max(mag(slaveTractionForMaster_));
    }
    reduce(maxMagMasterTraction, maxOp<scalar>());
    reduce(numSlipFaces, sumOp<int>());
    reduce(numStickFaces, sumOp<int>());

    // master writes to contact info file
    if (Pstream::master() && (contactIterNum_ %  infoFreq_ == 0))
    {
        OFstream& contactFile = *contactFilePtr_;
        int width = 20;
        contactFile
            << mesh.time().value();
        contactFile.width(width);
        contactFile
            << contactIterNum_;
        contactFile.width(width);
        contactFile.width(width);
        contactFile
            << numSlipFaces;
        contactFile.width(width);
        contactFile
            << numStickFaces;
        contactFile.width(width);
        contactFile
            << maxMagMasterTraction << endl;
    }
}


void Foam::dirichletNeumannFriction::writeDict(Ostream& os) const
{
    word keyword(name() + "FrictionModelDict");
    os.writeKeyword(keyword)
        << frictionContactModelDict_;
}


// ************************************************************************* //
