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
#include "PatchToPatchInterpolationTemplate.H"
#include "ggiInterpolation.H"
#include "tractionBoundaryGradient.H"
#include "fvc.H"
#include "primitivePatchInterpolation.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dirichletNeumannFriction, 0);
    addToRunTimeSelectionTable
    (frictionContactModel, dirichletNeumannFriction, dictionary);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


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
  oldSlip_(slaveDisp_.size(), vector::zero),
  slaveTraction_(mesh().boundaryMesh()[slavePatchID].size(), vector::zero),
  oldSlaveTraction_(slaveTraction_),
  slaveValueFrac_(mesh_.boundaryMesh()[slavePatchID].size(), symmTensor::zero),
  oldSlaveValueFrac_(slaveValueFrac_),
  relaxationFactor_
  (readScalar(frictionContactModelDict_.lookup("relaxationFactor"))),
  contactIterNum_(0),
  infoFreq_(readInt(frictionContactModelDict_.lookup("infoFrequency"))),
  oscillationCorr_(frictionContactModelDict_.lookup("oscillationCorrection")),
  smoothingSteps_(readInt(frictionContactModelDict_.lookup("smoothingSteps"))),
  oldStickSlip_(slaveDisp_.size(), 0.0),
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
      contactFile << "relaxationFactor";
      contactFile.width(width);
      contactFile << "slipFaces";
      contactFile.width(width);
      contactFile << "stickFaces";
      contactFile.width(width);
      contactFile << "maxMagSlaveTraction" << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

  void Foam::dirichletNeumannFriction::correct
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
    // to get local masterDUInterpToSlave (i.e. masterDU interpolated
    // to the slave)
    // so the method is:
    // create global masterDU field
    // interpolate global masterDU from master global face zone to
    // slave global zone
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
    FatalError << "dirichletNeumannFriction::correct()\n"
      " The displacement field must be called U or DU"
           << exit(FatalError);
      }

    // put local masterDU into globalMasterDU
    const label masterPatchStart
      = mesh.boundaryMesh()[masterPatchIndex].start();
    forAll(masterDU, i)
      {
          globalMasterDU[
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
    FatalError << "dirichletNeumannFriction::correct()\n"
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

    // Now masterDUInterpToSlave should have masterDU interpolated to the slave

    // Calculate current slave shear traction from the normal gradient field
    const fvPatch& slavePatch = mesh.boundary()[slavePatchIndex];
    const fvPatchField<tensor>& gradField =
        slavePatch.lookupPatchField<volTensorField, tensor>
        ("grad(" + fieldName + ")");

    bool incremental(fieldName == "DU");

    vectorField slaveShearTraction
    (
        (I - sqr(slaveFaceNormals))
      & tractionBoundaryGradient::traction
        (
            gradField,
            fieldName,
            "U",
            slavePatch,
            orthotropic,
            nonLinearGeometry::nonLinearNames_[nonLinear],
            incremental
        )
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
    const vectorField prevSlaveShearDisp =
        (I - sqr(slaveFaceNormals))
        & dispField.boundaryField()[slavePatchIndex];
    label numSlipFaces = 0;
    label numStickFaces = 0;
    scalarField& stickSlip = stickSlipFaces();
    const scalarField magSlavePressure = mag(slavePressure);
    scalar maxMagSlavePressure = 0.0;
    if (slavePressure.size() > 0) maxMagSlavePressure = max(magSlavePressure);
    reduce(maxMagSlavePressure, maxOp<scalar>());

    // slip is the difference between the master tangential DU
    // and slave tangential DU
    vectorField slip =
      (I - sqr(slaveFaceNormals)) &
      ( slaveDU - masterDUInterpToSlave);

    // under-relax the slip
    slip = relaxationFactor_*slip + (1.0 - relaxationFactor_)*oldSlip_;
    oldSlip_ = slip;

    vectorField slipDir = slip/(mag(slip)+SMALL);

    const scalar pressureTol = 1e-3*maxMagSlavePressure;


    // first we calculate slaveDisp assuming every face is sticking
    // and we calculate slaveTraction assuming every face is sliding
    // then we use the slaveValueFrac to set the face to slip/stick
    // depending on the slip function
    slaveDisp_ = -slip + prevSlaveShearDisp;
    //slaveDisp_ = -slip + slaveDisp_; // see if convergence is better

    scalarField slipTrac = frictionLawPtr_->slipTraction(magSlavePressure);

    // new
    forAll(slaveDisp_, facei)
      {
    if (mag(slavePressure[facei]) < pressureTol)
      {
        // not in contact
            slaveDisp_[facei] = prevSlaveShearDisp[facei];
            slaveTraction_[facei] = vector::zero;
            slaveValueFrac_[facei] = symmTensor::zero;

            stickSlip[facei] = -1;
      }
    else if (
        (mag(slaveShearTraction[facei]) > 0.999*slipTrac[facei])
        &&
        ((slip[facei] & slaveShearTraction[facei]) < 0.0) // opposite directions
        )
       {
           // slip
           slaveDisp_[facei] =
               -slip[facei] + prevSlaveShearDisp[facei]; // better convergence
         slaveTraction_[facei] = -slipDir[facei] * slipTrac[facei];
         slaveValueFrac_[facei] = symmTensor::zero;
             numSlipFaces++;

         stickSlip[facei] = 0;
       }
    else
       {
         // stick
             slaveDisp_[facei] = -slip[facei] + prevSlaveShearDisp[facei];
         slaveTraction_[facei] = -slipDir[facei] * slipTrac[facei];
         //slaveTraction_[facei] = slaveShearTraction[facei];
         //slaveTraction_[facei] = vector::zero;
         slaveValueFrac_[facei] = (I - sqr(slaveFaceNormals[facei]));
             numStickFaces++;

         stickSlip[facei] = 2;
       }
      }

    // slaveTraction_ = -slipDir * slipTrac;

    // // slipFunc is 1.0 for slip and 0.0 for stick
    // // scalarField slipFunc = mag(slaveShearTraction) - 0.999*slipTrac;
    // scalarField slipFunc = mag(slaveShearTraction) - slipTrac;
    // //slipFunc /= mag(slipFunc+SMALL); // bugfix add SMALL
    // //slipFunc = max(slipFunc, 0.0);
    // //const scalar slipStressTol = 1e-4*gMax(mag(slaveShearTraction));
    // const scalar pressureTol = 1e-3*maxMagSlavePressure;
    // forAll(slipFunc, facei)
    //   {
    //      if (slipFunc[facei] > pressureTol)
    //        {
    //          slipFunc[facei] = 1.0;
    //        }
    //      else
    //        {
    //          slipFunc[facei] = 0.0;
    //        }
    //   }

    // // if the slip and trac are in the same direction then we will change
    // // the face from slip to stick
    // // {
    // //   scalarField changeFace = slip & slaveTraction_;
    // //   forAll(changeFace, facei)
    // //   {
    // //     if (changeFace[facei] > SMALL)
    // //       {
    // //         slipFunc[facei] = 0.0;
    // //       }
    // //   }
    // // }

    // // stickSlip
    // // -1 for faces not in contact
    // // 0 for slipping faces
    // // 1 for stick faces
    // stickSlip = (1.0 - slipFunc);
    // // const scalar pressureTol = 1e-3*maxMagSlavePressure;
    // forAll(slip, facei)
    //   {
    //      if (magSlavePressure[facei] < pressureTol)
    //        {
    //          stickSlip[facei] = -1;
    //        }
    //   }

    // //slaveValueFrac_ = (1.0 - slipFunc) * (I - sqr(slaveFaceNormals));
    // slaveValueFrac_ = max(stickSlip, 0.0) * (I - sqr(slaveFaceNormals));

    // // set valueFace to zero for faces not in contact
    // // and also reset traction on sticking faces and displacement
    // // on slipping faces
    // forAll(slip, facei)
    //   {
    //      if (magSlavePressure[facei] < pressureTol)
    //        {
    //          slaveDisp_[facei] = prevSlaveShearDisp[facei];
    //          slaveTraction_[facei] = vector::zero;
    //          slaveValueFrac_[facei] = symmTensor::zero;
    //          if ( mag(stickSlip[facei] + 1) > SMALL)
    //            Info << "face " << facei << "changed to tracFree" << endl;

    //          stickSlip[facei] = -1;
    //        }
    //      else if (slipFunc[facei] > SMALL)
    //        {
    //          slaveDisp_[facei] = prevSlaveShearDisp[facei];
    //          numSlipFaces++;

    //          if ( mag(stickSlip[facei]) > SMALL)
    //            Info << "face " << facei << "changed to slip" << endl;

    //          stickSlip[facei] = 0;
    //        }
    //      else
    //        {
    //          //slaveTraction_[facei] = slaveShearTraction[facei];
    //          numStickFaces++;

    //          if ( mag(stickSlip[facei] - 1) > SMALL)
    //            Info << "face " << facei << "changed to stick" << endl;

    //          stickSlip[facei] = 1;
    //        }
    //   }

    // correct oscillations
    if (oscillationCorr_)
      {
    //correctOscillations(slaveFaceZonePatch);

    // interpolate face values to points then interpolate back
    // this essentially smooths the field
    primitivePatchInterpolation localSlaveInterpolator
        (mesh.boundaryMesh()[slavePatchIndex]);
    vectorField slaveDispPoints
        (mesh.boundaryMesh()[slavePatchIndex].nPoints(), vector::zero);

    for (int i=0; i<smoothingSteps_; i++)
      {
        slaveDispPoints =
            localSlaveInterpolator.faceToPointInterpolate<vector>(slaveDisp_);
        slaveDisp_ =
            localSlaveInterpolator.pointToFaceInterpolate<vector>
            (slaveDispPoints);

        // make sure no normal component
        slaveDisp_ = (I - sqr(slaveFaceNormals)) & slaveDisp_;
      }
      }

    // under-relax traction
    slaveDisp_ =
        relaxationFactor_*slaveDisp_
        + (1.0 - relaxationFactor_)*prevSlaveShearDisp; //oldSlaveDisp_;
    oldSlaveDisp_ = slaveDisp_;
    slaveValueFrac_ =
        relaxationFactor_*slaveValueFrac_
        + (1.0 - relaxationFactor_)*oldSlaveValueFrac_;
    oldSlaveValueFrac_ = slaveValueFrac_;
    slaveTraction_ =
        relaxationFactor_*slaveTraction_
        + (1.0 - relaxationFactor_)*oldSlaveTraction_;
    oldSlaveTraction_ = slaveTraction_;
    stickSlip =
        relaxationFactor_*stickSlip + (1.0 - relaxationFactor_)*oldStickSlip_;
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
        maxMagMasterTraction = max(mag(slaveTraction_));
    }
    reduce(maxMagMasterTraction, maxOp<scalar>());
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
    contactFile << relaxationFactor_;
    contactFile.width(width);
    contactFile << numSlipFaces;
    contactFile.width(width);
    contactFile << numStickFaces;
    contactFile.width(width);
    contactFile << maxMagMasterTraction << endl;
      }
  }



  void Foam::dirichletNeumannFriction::writeDict(Ostream& os) const
  {
    word keyword(name()+"FrictionModelDict");
    os.writeKeyword(keyword)
      << frictionContactModelDict_;
  }


// ************************************************************************* //
