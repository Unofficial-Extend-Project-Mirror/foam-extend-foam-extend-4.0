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

#include "standardPenaltyFriction.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include "patchToPatchInterpolation.H"
#include "ggiInterpolation.H"
#include "constitutiveModel.H"
#include "fvc.H"
#include "primitivePatchInterpolation.H"
#include "solidContactFvPatchVectorField.H"
//#include "solidGeneralContactFvPatchVectorField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(standardPenaltyFriction, 0);
    addToRunTimeSelectionTable
    (
        frictionContactModel,
        standardPenaltyFriction,
        dictionary
    );
}


// * * * * * * * * * * * Private Members Functions * * * * * * * * * * * * * //

void Foam::standardPenaltyFriction::calcFrictionPenaltyFactor()
{
    // set penalty factor using a similar method to the normal
    // contact where we approx penaltyFactor from mechanical properties
    // this can then be scaled using the penaltyScale
    const label masterPatchIndex =  masterPatchID();
    const label slavePatchIndex =  slavePatchID();
    const constitutiveModel& mechanical =
        mesh_.objectRegistry::lookupObject<constitutiveModel>
        (
            "rheologyProperties"
        );
    scalarField masterMu =
        mechanical.mu()().boundaryField()[masterPatchIndex];
    scalarField slaveMu =
        mechanical.mu()().boundaryField()[slavePatchIndex];

    // avarage contact patch shear modulus
    scalar shearModulus = 0.5*(gAverage(masterMu) + gAverage(slaveMu));

    // average contact patch face area
    scalar masterMagSf =
        gAverage(mesh_.magSf().boundaryField()[masterPatchIndex]);
    scalar slaveMagSf =
        gAverage(mesh_.magSf().boundaryField()[slavePatchIndex]);
    scalar faceArea = 0.5*(masterMagSf + slaveMagSf);

    // average contact patch cell volume
    scalarField masterV(mesh_.boundary()[masterPatchIndex].size(), 0.0);
    scalarField slaveV(mesh_.boundary()[slavePatchIndex].size(), 0.0);
    const volScalarField::DimensionedInternalField& V = mesh_.V();
    {
        const unallocLabelList& faceCells =
            mesh_.boundary()[masterPatchIndex].faceCells();
        forAll(mesh_.boundary()[masterPatchIndex], facei)
        {
            masterV[facei] = V[faceCells[facei]];
        }
    }
    {
        const unallocLabelList& faceCells =
            mesh_.boundary()[slavePatchIndex].faceCells();
        forAll(mesh_.boundary()[slavePatchIndex], facei)
        {
            slaveV[facei] = V[faceCells[facei]];
        }
    }
    scalar cellVolume = 0.5*(gAverage(masterV) + gAverage(slaveV));

    // approximate penalty factor based on Hallquist et al.
    // we approximate penalty factor for traction instead of force
    frictionPenaltyFactorPtr_ =
        new scalar(frictionPenaltyScale_*shearModulus*faceArea/cellVolume);

    Info<< "    friction penalty factor: " << *frictionPenaltyFactorPtr_
        << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::standardPenaltyFriction::standardPenaltyFriction
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
    slaveTraction_(mesh().boundaryMesh()[slavePatchID].size(), vector::zero),
    oldSlaveTraction_(mesh().boundaryMesh()[slavePatchID].size(), vector::zero),
    slaveValueFrac_
    (
        mesh().boundaryMesh()[slavePatchID].size(), symmTensor::zero
    ),
    slip_(slaveDisp_.size(), vector::zero),
    frictionPenaltyFactorPtr_(NULL),
    frictionPenaltyScale_
    (
        readScalar(frictionContactModelDict_.lookup("penaltyScale"))
    ),
    relaxFac_(readScalar(frictionContactModelDict_.lookup("relaxationFactor"))),
    contactIterNum_(0),
    infoFreq_(readInt(frictionContactModelDict_.lookup("infoFrequency"))),
    oscillationCorr_(frictionContactModelDict_.lookup("oscillationCorrection")),
    smoothingSteps_
    (
        readInt(frictionContactModelDict_.lookup("smoothingSteps"))
    ),
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

    // master proc open contact info file
    if (Pstream::master())
    {
        word masterName = mesh_.boundary()[masterPatchID].name();
        word slaveName = mesh_.boundary()[slavePatchID].name();
        fileName contactFileDir = "contact";
        mkDir(contactFileDir);
        contactFilePtr_ =
            new OFstream
            (
                contactFileDir/
                "frictionContact_"+masterName+"_"+slaveName+".txt"
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
            << "penaltyScale";
        contactFile.width(width);
        contactFile
            << "slipFaces";
        contactFile.width(width);
        contactFile
            << "stickFaces";
        contactFile.width(width);
        contactFile
            << "maxMagSlaveTraction";
        contactFile.width(width);
        contactFile
            << "maxMagSlip" << endl;
    }
}


// Construct as a copy
Foam::standardPenaltyFriction::standardPenaltyFriction
(
    const standardPenaltyFriction& fm
)
:
    frictionContactModel(fm),
    frictionContactModelDict_(fm.frictionContactModelDict_),
    frictionLawPtr_(fm.frictionLawPtr_->clone().ptr()),
    mesh_(fm.mesh_),
    slaveDisp_(fm.slaveDisp_),
    slaveTraction_(fm.slaveTraction_),
    oldSlaveTraction_(fm.oldSlaveTraction_),
    slaveValueFrac_(fm.slaveValueFrac_),
    slip_(fm.slip_),
    frictionPenaltyFactorPtr_(NULL),
    frictionPenaltyScale_(fm.frictionPenaltyScale_),
    relaxFac_(fm.relaxFac_),
    contactIterNum_(fm.contactIterNum_),
    infoFreq_(fm.infoFreq_),
    oscillationCorr_(fm.oscillationCorr_),
    smoothingSteps_(fm.smoothingSteps_),
    contactFilePtr_(NULL)
{
    if (fm.frictionPenaltyFactorPtr_)
    {
        frictionPenaltyFactorPtr_ = new scalar(*fm.frictionPenaltyFactorPtr_);
    }

    if (fm.contactFilePtr_)
    {
        contactFilePtr_ = new OFstream(*fm.contactFilePtr_);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::standardPenaltyFriction::correct
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
    // slave); so the method is:
    // create global masterDU field;
    // interpolate global masterDU from master global face zone to slave global
    // zone;
    // then find local masterDUInterpToSlave from the global interpolated field.

    vectorField masterDUInterpToSlave
        (
            mesh.boundaryMesh()[slavePatchIndex].size(), vector::zero
        );

    // global master DU
    vectorField globalMasterDU(masterFaceZonePatch.size(), vector::zero);

    // lookup current displacement field
    const volVectorField& dispField =
        mesh.objectRegistry::lookupObject<volVectorField>(fieldName);

    // local master and slave DU increment
    vectorField masterDU = dispField.boundaryField()[masterPatchIndex];
    vectorField slaveDU = dispField.boundaryField()[slavePatchIndex];

    //if (fieldName == "U")
    if (fieldName == "D")
    {
        // lookup old U
        // const volVectorField& dispOldField =
        //    mesh.objectRegistry::lookupObject<volVectorField>(fieldName+"_0");

        // // subtract old U
        // masterDU -= dispOldField.boundaryField()[masterPatchIndex];
        // slaveDU -= dispOldField.boundaryField()[slavePatchIndex];
        // philipc fix: use oldTime function
        masterDU -= dispField.oldTime().boundaryField()[masterPatchIndex];
        slaveDU -= dispField.oldTime().boundaryField()[slavePatchIndex];
    }
    //else if (fieldName != "DU")
    else if (fieldName != "DD")
    {
        FatalErrorIn("standardPenaltyFriction::correct()")
            << " The displacement field must be called U or DU"
            << abort(FatalError);
    }

    // put local masterDU into globalMasterDU
    const label masterPatchStart
        = mesh.boundaryMesh()[masterPatchIndex].start();
    forAll(masterDU, i)
    {
        globalMasterDU
        [
            mesh.faceZones()[masterFaceZoneID()].whichFace(masterPatchStart + i)
        ] = masterDU[i];
    }

    //- exchange parallel data
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
        FatalErrorIn("standardPenaltyFriction::correct()")
            << "interpolationMethod " << interpolationMethod << " not known\n"
            << "interpolationMethod must be patchToPatch or ggi"
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
                mesh.faceZones()[
                    slaveFaceZoneID()
                    ].whichFace(slavePatchStart + i)
            ];
    }

    // Now masterDUInterpToSlave should be masterDU interpolated to the slave

    // Calculate slave shear traction increments
    const scalarField magSlavePressure = mag(slavePressure);
    label numSlipFaces = 0;
    label numStickFaces = 0;
    scalarField& stickSlip = stickSlipFaces();
    const scalarField oldStickSlip = stickSlip;
    scalar frictionPenaltyFac = frictionPenaltyFactor();
    scalar maxMagSlip = 0.0;
    scalarField slipTraction(magSlavePressure.size(), 0.0);

    // lookup areaInContact
    const scalarField* areaInContactPtr = NULL;

    if
    (
        isA<solidContactFvPatchVectorField>
        (
            dispField.boundaryField()[masterPatchIndex]
        )
    )
    {
        const solidContactFvPatchVectorField& DUpatch =
            refCast<const solidContactFvPatchVectorField>
            (
                dispField.boundaryField()[masterPatchIndex]
            );
        areaInContactPtr = &(DUpatch.normalModel().areaInContact());
    }
    // else if
    // (
    //     isA<solidGeneralContactFvPatchVectorField>
    //     (
    //         dispField.boundaryField()[masterPatchIndex]
    //     )
    // )
    // {
    //     const solidGeneralContactFvPatchVectorField& DUpatch =
    //         refCast<const solidGeneralContactFvPatchVectorField>
    //         (
    //             //dispField.boundaryField()[masterPatchIndex]
    //             dispField.boundaryField()[slavePatchIndex]
    //         );

    //     const label shadowI = DUpatch.findShadowID(masterPatchIndex);

    //     areaInContactPtr = &(DUpatch.normalModel(shadowI).areaInContact());
    // }
    else
    {
        FatalErrorIn("standardPenaltyFriction::correct()")
            << "The contact patch must of type solidContact or "
            << "solidGeneralContact" << abort(FatalError);
    }

    const scalarField& areaInContact = *areaInContactPtr;


    forAll(magSlavePressure, faceI)
    {
        // there can only be a frictional tangential force when there is
        // a positive pressure
        // if (magSlavePressure[faceI] > pressureTol) //SMALL)
        if (areaInContact[faceI] > SMALL)
        {
            // Compute slip
            //- we need the difference of DU between the master and slave
            slip_[faceI] = slaveDU[faceI] - masterDUInterpToSlave[faceI];
            //- the shear traction direction is got by removing the normal
            // component of the DU
            //- (I - sqr(n)) removes the normal
            //- sqr(n) would remove the shear
            slip_[faceI] = (I - sqr(slaveFaceNormals[faceI])) & slip_[faceI];

            slaveTraction_[faceI] = -frictionPenaltyFac*slip_[faceI];
            const scalar magSlip = mag(slip_[faceI]);
            maxMagSlip = max(maxMagSlip, magSlip);

            // traction to cause slipping
            slipTraction[faceI] =
                frictionLawPtr_->slipTraction(magSlavePressure[faceI]);

            scalar slipFunc = mag(slaveTraction_[faceI]) - slipTraction[faceI];
            if (slipFunc > SMALL)
            {
                // analogous to plasticity
                // slip is a combination of elastic slip
                // and plastic slip
                // elastic slip should be zero but is finite due to
                // penalty stiffness
                // plastic slip is the permanent deformation
                slaveTraction_[faceI] =
                    slipTraction[faceI]*(-slip_[faceI]/magSlip);

                numSlipFaces++;
                stickSlip[faceI] = 1;
            }
            else
            {
                numStickFaces++;
                stickSlip[faceI] = 2;
            }
        }
        // no friction if pressure is negative or zero or face is not in contact
        else
        {
            slaveTraction_[faceI] = vector::zero;
            slipTraction[faceI] = 0.0;
            stickSlip[faceI] = 0;
        }
    }

    // stickSlip field is just for visualisation but we will under-relax it to
    // allow us to see if a face is jumping between stick and slip
    stickSlip = relaxFac_*stickSlip + (1.0 - relaxFac_)*oldStickSlip;

    if (oscillationCorr_)
    {
        // interpolate face values to points then interpolate back
        // this essentially smooths the field
        primitivePatchInterpolation localSlaveInterpolator
            (
                mesh.boundaryMesh()[slavePatchIndex]
            );
        vectorField slaveTracPoints
            (
                mesh.boundaryMesh()[slavePatchIndex].nPoints(), vector::zero
            );

        for (int i=0; i<smoothingSteps_; i++)
        {
            slaveTracPoints =
                localSlaveInterpolator.faceToPointInterpolate<vector>
                (
                    slaveTraction_
                );
            slaveTraction_ =
                localSlaveInterpolator.pointToFaceInterpolate<vector>
                (
                    slaveTracPoints
                );
        }
    }

    // Under-relax traction
    slaveTraction_ =
        relaxFac_*slaveTraction_ + (1.0 - relaxFac_)*oldSlaveTraction_;
    oldSlaveTraction_ = slaveTraction_;

    scalar maxMagSlaveTraction = 0.0;
    if (slaveTraction_.size() > 0)
    {
        maxMagSlaveTraction = max(mag(slaveTraction_));
    }
    reduce(maxMagSlaveTraction, maxOp<scalar>());
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
        contactFile
            << frictionPenaltyScale_;
        contactFile.width(width);
        contactFile
            << numSlipFaces;
        contactFile.width(width);
        contactFile
            << numStickFaces;
        contactFile.width(width);
        contactFile
            << maxMagSlaveTraction;
        contactFile.width(width);
        contactFile
            << maxMagSlip << endl;
    }
}


void Foam::standardPenaltyFriction::writeDict(Ostream& os) const
{
    word keyword(name()+"FrictionModelDict");
    os.writeKeyword(keyword)
        << frictionContactModelDict_;
}


// ************************************************************************* //
