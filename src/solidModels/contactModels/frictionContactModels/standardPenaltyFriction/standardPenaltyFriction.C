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
#include "PatchToPatchInterpolationTemplate.H"
#include "ggiInterpolation.H"
#include "constitutiveModel.H"
#include "fvc.H"
#include "primitivePatchInterpolation.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(standardPenaltyFriction, 0);
    addToRunTimeSelectionTable
    (frictionContactModel, standardPenaltyFriction, dictionary);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


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
  slaveValueFrac_(mesh().boundaryMesh()[slavePatchID].size(), symmTensor::zero),
  frictionPenaltyFactorPtr_(NULL),
  frictionPenaltyScale_
  (readScalar(frictionContactModelDict_.lookup("penaltyScale"))),
  relaxFac_(readScalar(frictionContactModelDict_.lookup("relaxationFactor"))),
  contactIterNum_(0),
  infoFreq_(readInt(frictionContactModelDict_.lookup("infoFrequency"))),
  oscillationCorr_(frictionContactModelDict_.lookup("oscillationCorrection")),
  smoothingSteps_(readInt(frictionContactModelDict_.lookup("smoothingSteps"))),
  contactFilePtr_(NULL)
{
  // create friction law
  frictionLawPtr_ = frictionLaw::New(
      frictionContactModelDict_.lookup("frictionLaw"),
      frictionContactModelDict_
      ).ptr();

  // master proc open contact info file
  if (Pstream::master())
    {
      word masterName = mesh_.boundary()[masterPatchID].name();
      word slaveName = mesh_.boundary()[slavePatchID].name();
      contactFilePtr_ =
          new OFstream
          (fileName("frictionContact_"+masterName+"_"+slaveName+".txt"));
      OFstream& contactFile = *contactFilePtr_;
      int width = 20;
      contactFile << "time";
      contactFile.width(width);
      contactFile << "iterNum";
      contactFile.width(width);
      contactFile << "penaltyScale";
      contactFile.width(width);
      contactFile << "slipFaces";
      contactFile.width(width);
      contactFile << "stickFaces";
      contactFile.width(width);
      contactFile << "maxMagSlaveTraction" << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

  void Foam::standardPenaltyFriction::correct
  (
   const vectorField& slavePressure,
   const PrimitivePatch<face, List, pointField>& masterFaceZonePatch,
   const PrimitivePatch<face, List, pointField>& slaveFaceZonePatch,
   const intersection::algorithm alg,
   const intersection::direction dir,
   const word interpolationMethod,
   const word fieldName,
   const Switch orthotropic,
   const word nonLinear,
   const vectorField& slaveFaceNormals
   )
  {
    const fvMesh& mesh = mesh_;
    const label slavePatchIndex = slavePatchID();
    const label masterPatchIndex = masterPatchID();
    contactIterNum_++;

    // we have local masterDU and we want to interpolate it to the slave
    // to get local masterDUInterpToSlave (i.e. masterDU interpolated to
    // the slave)
    // so the method is:
    // create global masterDU field
    // interpolate global masterDU from master global face zone to slave
    // global zone
    // then find local masterDUInterpToSlave from the global interpolated field

    vectorField masterDUInterpToSlave
        (mesh.boundaryMesh()[slavePatchIndex].size(), vector::zero);

    // global master DU
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
      mesh.objectRegistry::lookupObject<volVectorField>(fieldName+"_0");

    // subtract old U
    masterDU -= dispOldField.boundaryField()[masterPatchIndex];
    slaveDU -= dispOldField.boundaryField()[slavePatchIndex];
      }
    else if (fieldName != "DU")
      {
    FatalError << "iterativePenaltyFunction::correct()\n"
      " The displacement field must be called U or DU"
           << exit(FatalError);
      }

    // put local masterDU into globalMasterDU
    const label masterPatchStart
      = mesh.boundaryMesh()[masterPatchIndex].start();
    forAll(masterDU, i)
      {
          globalMasterDU
              [
                  mesh.faceZones()[masterFaceZoneID()
                      ].whichFace(masterPatchStart + i)] =
              masterDU[i];
      }
    //- exchange parallel data
    // sum because each face is only on one proc
    reduce(globalMasterDU, sumOp<vectorField>());

    // globalMasterDU is interpolated to the slave
    vectorField globalMasterDUInterpToSlave
        (slaveFaceZonePatch.size(), vector::zero);

    // interpolate DU from master to slave using inverse distance or ggi
    if (interpolationMethod == "patchToPatch")
      {
    PatchToPatchInterpolation<
      PrimitivePatch<
          face, List, pointField
          >, PrimitivePatch<face, List, pointField>
      > masterToSlavePatchToPatchInterpolator
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
    GGIInterpolation<
      PrimitivePatch<
          face, List, pointField
          >, PrimitivePatch< face, List, pointField >
      > masterToSlaveGgiInterpolator
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
    FatalError << "iterativePenaltyFunction::correct()\n"
      "interpolationMethod " << interpolationMethod << " not known\n"
      "interpolationMethod must be patchToPatch or ggi"
           << exit(FatalError);
      }

    // now put global back into local
    const label slavePatchStart
      = mesh.boundaryMesh()[slavePatchIndex].start();

    forAll(masterDUInterpToSlave, i)
      {
    masterDUInterpToSlave[i] =
      globalMasterDUInterpToSlave
      [
       mesh.faceZones()[slaveFaceZoneID()].whichFace(slavePatchStart + i)
       ];
      }

    // Now masterDUInterpToSlave should be masterDU interpolated to the slave

    // calculate slave shear traction increments
    const scalarField magSlavePressure = mag(slavePressure);
    label numSlipFaces = 0;
    label numStickFaces = 0;
    scalarField& stickSlip = stickSlipFaces();
    scalar frictionPenaltyFac = frictionPenaltyFactor();
    scalar maxMagSlip = 0.0; // of stick faces
    //scalar slipTol = 5e-8;
    scalar pressureTol = 1e-3*gMax(magSlavePressure);
    scalarField slipTraction(magSlavePressure.size(), 0.0);

    forAll(magSlavePressure, faceI)
      {
    // there can only be a frictional tangential force when there is
    // a positive pressure
    if (magSlavePressure[faceI] > pressureTol) //SMALL)
      {
        // compute slip
        //- we need the difference of DU between the master and slave
        vector slip =
          slaveDU[faceI] - masterDUInterpToSlave[faceI];
        //- the shear traction direction is got by removing the normal
        // component of the DU
        //- (I - sqr(n)) removes the normal
        //- sqr(n) would remove the shear
        slip = (I - sqr(slaveFaceNormals[faceI])) & slip;

        //slaveTrac = -frictionPenaltyFac * slip;
        slaveTraction_[faceI] = -frictionPenaltyFac * slip;
        maxMagSlip = max(maxMagSlip, mag(slip));

        // traction to cause slipping
        slipTraction[faceI] =
            frictionLawPtr_->slipTraction(magSlavePressure[faceI]);
        scalar slipFunc = mag(slaveTraction_[faceI]) - 0.99*slipTraction[faceI];
        if (slipFunc > 0)
          {
        // analogous to plasticity
        // slip is a combination of elastic slip
        // and plastic slip
        // elastic slip should be zero but is finite due to
        // penalty stiffness
        // plastic slip is the permanent deformation
        slaveTraction_[faceI] = slipTraction[faceI]*(-slip/mag(slip));

        numSlipFaces++;
        stickSlip[faceI] = 1;
          }
        else
          {
        numStickFaces++;
        stickSlip[faceI] = 2;
          }
      }
    // no friction if pressure is negative or zero
    else
      {
        slaveTraction_[faceI] = vector::zero;
        slipTraction[faceI] = 0.0;
        stickSlip[faceI] = 0;
      }
      }

    if (oscillationCorr_)
      {
        //correctOscillations(slaveTraction_, slipTraction, slaveFaceZonePatch);

    // interpolate face values to points then interpolate back
    // this essentially smooths the field
    primitivePatchInterpolation localSlaveInterpolator
        (mesh.boundaryMesh()[slavePatchIndex]);
    vectorField slaveTracPoints
        (mesh.boundaryMesh()[slavePatchIndex].nPoints(), vector::zero);

    for (int i=0; i<smoothingSteps_; i++)
      {
        slaveTracPoints =
            localSlaveInterpolator.faceToPointInterpolate<vector>
            (slaveTraction_);
        slaveTraction_ =
            localSlaveInterpolator.pointToFaceInterpolate<vector>
            (slaveTracPoints);
      }
      }

    // under-relax traction
    slaveTraction_ =
        relaxFac_*slaveTraction_ + (1.0 - relaxFac_)*oldSlaveTraction_;

    // apply extra relaxation to faces on the border of stick slip zones
    // - only serial - testing
    // forAll(slaveTraction_, facei)
    //   {
    //      scalar relax = false;
    //      forAll(slaveFaceZonePatch.faceFaces()[facei], ffi)
    //        {
    //          label faceFace = slaveFaceZonePatch.faceFaces()[facei][ffi];

    //          if (mag(stickSlip[facei] - stickSlip[faceFace]) > SMALL
    //             &&
    //             mag(stickSlip[facei] - 1.0) < SMALL)
    //            {
    //          relax = true;
    //          break;
    //            }
    //        }
    //      if (relax)
    //        {
    //          slaveTraction_[facei] =
    //            relaxFac_*slaveTraction_[facei]
    // + (1.0 - relaxFac_)*oldSlaveTraction_[facei];
    //        }
    //   }

    oldSlaveTraction_ = slaveTraction_;

    scalar maxMagSlaveTraction = gMax(mag(slaveTraction_));
    reduce(numSlipFaces, sumOp<int>());
    reduce(numStickFaces, sumOp<int>());

    // master writes to contact info file
    if (Pstream::master() && (contactIterNum_ %  infoFreq_ == 0))
      {
    OFstream& contactFile = *contactFilePtr_;
    int width = 20;
    contactFile << mesh.time().value();
    contactFile.width(width);
    contactFile << contactIterNum_;
    contactFile.width(width);
    contactFile << frictionPenaltyScale_;
    contactFile.width(width);
    contactFile << numSlipFaces;
    contactFile.width(width);
    contactFile << numStickFaces;
    contactFile.width(width);
    contactFile << maxMagSlaveTraction;
    contactFile.width(width);
    contactFile << maxMagSlip;
    contactFile << endl;
      }
  }



  void Foam::standardPenaltyFriction::calcFrictionPenaltyFactor()
  {
    // set penalty factor using a similar method to the normal
    // contact where we approx penaltyFactor from mechanical properties
    // this can then be scaled using the penaltyScale
    const label masterPatchIndex =  masterPatchID();
    const label slavePatchIndex =  slavePatchID();
    const constitutiveModel& rheology =
      mesh_.objectRegistry::lookupObject<constitutiveModel>
        ("rheologyProperties");
    scalarField masterMu =
      rheology.mu()().boundaryField()[masterPatchIndex];
    scalarField slaveMu =
      rheology.mu()().boundaryField()[slavePatchIndex];

    // avarage contact patch shear modulus
    scalar shearModulus = 0.5*(gAverage(masterMu)+gAverage(slaveMu));

    // average contact patch face area
    scalar masterMagSf =
        gAverage(mesh_.magSf().boundaryField()[masterPatchIndex]);
    scalar slaveMagSf =
        gAverage(mesh_.magSf().boundaryField()[slavePatchIndex]);
    scalar faceArea = 0.5*(masterMagSf+slaveMagSf);

    // average contact patch cell volume
    scalarField masterV(mesh_.boundary()[masterPatchIndex].size(), 0.0);
    scalarField slaveV(mesh_.boundary()[slavePatchIndex].size(), 0.0);
    const volScalarField::DimensionedInternalField & V = mesh_.V();
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
    scalar cellVolume = 0.5*(gAverage(masterV)+gAverage(slaveV));

    // approximate penalty factor based on Hallquist et al.
    // we approximate penalty factor for traction instead of force
    frictionPenaltyFactorPtr_ =
        new scalar(frictionPenaltyScale_*shearModulus*faceArea/cellVolume);
  }

  void Foam::standardPenaltyFriction::writeDict(Ostream& os) const
  {
    word keyword(name()+"FrictionModelDict");
    os.writeKeyword(keyword)
      << frictionContactModelDict_;
  }


// ************************************************************************* //
