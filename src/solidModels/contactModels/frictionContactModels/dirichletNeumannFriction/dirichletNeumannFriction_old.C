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
  slaveTraction_(mesh().boundaryMesh()[slavePatchID].size(), vector::zero),
  slaveValueFrac_(mesh_.boundaryMesh()[slavePatchID].size(), symmTensor::zero),
  relaxationFactor_
  (readScalar(frictionContactModelDict_.lookup("relaxationFactor"))),
  contactIterNum_(0),
  infoFreq_(readInt(frictionContactModelDict_.lookup("infoFrequency"))),
  oscillationCorr_(frictionContactModelDict_.lookup("oscillationCorrection")),
  oscillationCorrFac_
  (readScalar(frictionContactModelDict_.lookup("oscillationCorrectionFactor"))),
  contactFilePtr_(NULL)
{
  // create friction law
  frictionLawPtr_ =
      frictionLaw::New(
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
    FatalError
        << "iterativePenaltyFunction::correct()\n"
        << "interpolationMethod " << interpolationMethod << " not known\n"
        << "interpolationMethod must be patchToPatch or ggi"
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

   // calculate current slave shear traction
    const fvPatch& slavePatch = mesh.boundary()[slavePatchIndex];
    const fvPatchField<tensor>& gradField =
      slavePatch.lookupPatchField<volTensorField, tensor>
        ("grad("+fieldName+")");
    vectorField slaveShearTraction =
      (I - sqr(slaveFaceNormals))
      &
      tractionBoundaryGradient().traction
      (
       gradField,
       fieldName,
       slavePatch,
       orthotropic,
       nonLinear
       );


    // calculate slave shear displacement increments
    const scalarField magSlavePressure = -1.0*(slaveFaceNormals&slavePressure);
    const volVectorField& oldSlaveDispField =
      mesh.objectRegistry::lookupObject<volVectorField>(fieldName);
    const vectorField& oldSlaveDisp =
        oldSlaveDispField.boundaryField()[slavePatchIndex];
    label numSlipFaces = 0;
    label numStickFaces = 0;
    scalarField& stickSlip = stickSlipFaces();

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
    const scalar maxMagSlavePressure = gMax(magSlavePressure);
    forAll(magSlavePressure, faceI)
      {
    // there can only be a frictional tangential force when there is
    // a positive pressure
    // if (magSlavePressure[faceI] > SMALL)
    if (magSlavePressure[faceI] > 1e-3*maxMagSlavePressure)
      {
        //scalar slipTrac = frictionCoeff_*magSlavePressure[faceI];
        scalar slipTrac =
            frictionLawPtr_->slipTraction(magSlavePressure[faceI]);

        // slipping faces
        if (mag(slaveShearTraction[faceI]) > (0.99*slipTrac) )
          {
        // direction of shear traction
        vector tracDir =
            slaveShearTraction[faceI] / mag(slaveShearTraction[faceI]);

        // slip is the difference between the master tangential DU
        // and slave tangential DU
        vector slip =
          (I - sqr(slaveFaceNormals[faceI])) &
          ( slaveDU[faceI] - masterDUInterpToSlave[faceI]);

        // if the slip and dir are in the same direction then
        // we will make this a
        // sticking face
        if ((tracDir & slip) > SMALL)
          {
            //Info << "face " << faceI << " flipping direction" << endl;
            numStickFaces++;
            stickSlip[faceI] = 2;

            // increment the slave shear displacement
            // we add an increment of shear disp to the slave faces
            // until there is no more slip
            slaveDisp_[faceI] =
              -1*relaxationFactor_ * slip;

            // slaveDisp_[faceI] is the correction to the disp so we
            // add on the original disp
            slaveDisp_[faceI] += oldSlaveDisp[faceI];
            // remove normal component
            slaveDisp_[faceI] =
                (I-sqr(slaveFaceNormals[faceI])) & slaveDisp_[faceI];

            // set slave valueFraction
            slaveValueFrac_[faceI] =
              relaxationFactor_*(I - sqr(slaveFaceNormals[faceI]))
              + (1.0 - relaxationFactor_)*slaveValueFrac_[faceI];

            // update traction as it is passed to the master
            slaveTraction_[faceI] =
              relaxationFactor_*slaveShearTraction[faceI]
              + (1-relaxationFactor_)*slaveTraction_[faceI];
          }
        // else we will limit the shear traction to slipTrac
        else
          {
            numSlipFaces++;
            stickSlip[faceI] = 1;

            // limit shear traction
            slaveTraction_[faceI] =
              relaxationFactor_*slipTrac*tracDir
              + (1-relaxationFactor_)*slaveTraction_[faceI];

            // update slave disp although it is not used for this face
            // while slipping
            slaveDisp_[faceI] =
                (I-sqr(slaveFaceNormals[faceI])) & oldSlaveDisp[faceI];

            // relax the slave valueFraction to zero
            //slaveValueFrac_[faceI] =
            (1.0 - relaxationFactor_)*slaveValueFrac_[faceI];
            slaveValueFrac_[faceI] = symmTensor::zero;
          }
          }
        // sticking faces
        else
          {
        numStickFaces++;
        stickSlip[faceI] = 2;

        // slip is the difference of the tangential DU between
        // the master and slave
        vector slip =
          (I - sqr(slaveFaceNormals[faceI])) &
          (slaveDU[faceI] - masterDUInterpToSlave[faceI]);

        // increment the slave shear displacement
        // we add an increment of shear disp to the slave faces
        // until there is no more slip
        slaveDisp_[faceI] = -1*relaxationFactor_*slip;
        // slaveDisp_[faceI] is the correction to the disp so we
        // add on the original disp
        slaveDisp_[faceI] += oldSlaveDisp[faceI];
        // remove normal component
        slaveDisp_[faceI] =
            (I-sqr(slaveFaceNormals[faceI])) & slaveDisp_[faceI];

        // set slave valueFraction
        slaveValueFrac_[faceI] =
          relaxationFactor_*(I - sqr(slaveFaceNormals[faceI]))
          + (1.0 - relaxationFactor_)*slaveValueFrac_[faceI];

        // update traction as it is passed to the master
        slaveTraction_[faceI] =
          relaxationFactor_*slaveShearTraction[faceI]
          + (1-relaxationFactor_)*slaveTraction_[faceI];
          }
      }
    // no friction if pressure is negative or zero
    else
      {
        stickSlip[faceI] = 0;
        // relax to zero
        slaveTraction_[faceI] = (1.0 - relaxationFactor_)*slaveTraction_[faceI];
        slaveValueFrac_[faceI] =
            (1.0 - relaxationFactor_)*slaveValueFrac_[faceI];
      }
      }

    // correct oscillations
    if (oscillationCorr_)
      {
    correctOscillations(slaveFaceZonePatch);
      }

    // get global values
    // in parallel, the log is poluted with warnings that
    // I am getting max of a list of size zero so
    // I will get the max of procs which have some
    // of the slave faces
    //scalar maxMagMasterTraction = gMax(mag(slaveTraction_))
    scalar maxMagMasterTraction = 0.0;
    if (slaveTraction_.size() > 0)
        maxMagMasterTraction = max(mag(slaveTraction_));
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

  void Foam::dirichletNeumannFriction::correctOscillations
  (
   const PrimitivePatch<face, List, pointField>& slaveFaceZonePatch
   )
  {
    // oscillations sometimes appear in contact shear displacements/tractions
    // so we will try to limit them here
    // we will weight the current face slaveDisp/Traction with
      // the average of the
    // neighbours using the weight oscillationCorrectionFactor_

    //Pout << "Applying contact shear oscillation correction..." << flush;

    const fvMesh& mesh = mesh_;
    const label slavePatchIndex = slavePatchID();
    const labelListList& faceFaces = slaveFaceZonePatch.faceFaces();
    const scalarField& stickSlip = stickSlipFaces();

    // create global fields
    //vectorField globalSlaveTraction(slaveFaceZonePatch.size(), vector::zero);
    vectorField globalSlaveDisp(slaveFaceZonePatch.size(), vector::zero);

    scalarField globalStickSlip(slaveFaceZonePatch.size(), 0.0);
    const label slavePatchStart
      = mesh.boundaryMesh()[slavePatchIndex].start();
    forAll(slaveTraction_, i)
      {
          globalSlaveDisp
              [
                  mesh.faceZones()[slaveFaceZoneID()
                      ].whichFace(slavePatchStart + i)] =
              slaveDisp_[i];
          globalStickSlip
              [
                  mesh.faceZones()[slaveFaceZoneID()
                      ].whichFace(slavePatchStart + i)] =
              stickSlip[i];
      }
    // sum because each face is only on one proc
    //reduce(globalSlaveTraction, sumOp<vectorField>());
    reduce(globalSlaveDisp, sumOp<vectorField>());
    reduce(globalStickSlip, sumOp<scalarField>());

    // smooth mag of slaveTraction with face face disps
    forAll(faceFaces, facei)
      {
    // only smooth sticking faces
    //if (mag(globalSlaveValueFrac[facei]) > SMALL)
    if (mag(globalStickSlip[facei] - 2.0) < SMALL)
      {
        //vector avTrac = vector::zero;
        vector avDisp = vector::zero;
        int numNei = 0;
        forAll(faceFaces[facei], ffi)
          {
        label faceFace = faceFaces[facei][ffi];

        // only include other sticking faces
        if ( mag(globalStickSlip[faceFace] - 2.0) < SMALL )
          {
            avDisp += globalSlaveDisp[faceFace];
            numNei++;
          }
          }

        // avTracMag /= numNei;
        //avTrac /= numNei;
        // if (numNei > 0)
        if (numNei > 1)
          {
        avDisp /= numNei;
          }
        else
          {
        avDisp = globalSlaveDisp[facei];
          }

        // if (numFaceFaces == 1)
        //   {
        //      // for corner/end faces, decrease the weight of the neighbours
        //      avTracMag += globalSlaveTraction[facei];
        //      avTracMag /= 2;
        //   }

        // weighted-average with face-faces
        globalSlaveDisp[facei] =
          oscillationCorrFac_*globalSlaveDisp[facei]
            + (1.0-oscillationCorrFac_)*avDisp;
      }
      }

    // convert global back to local
    forAll(slaveTraction_, facei)
      {
    // slaveTraction_[facei] =
    //   globalSlaveTraction
    //   [
    //    mesh.faceZones()[slaveFaceZoneID()].whichFace(slavePatchStart + facei)
    //    ];
    slaveDisp_[facei] =
      globalSlaveDisp
      [
       mesh.faceZones()[slaveFaceZoneID()].whichFace(slavePatchStart + facei)
       ];
      }

    //Pout << "\tdone" << endl;
  }


  void Foam::dirichletNeumannFriction::writeDict(Ostream& os) const
  {
    word keyword(name()+"FrictionModelDict");
    os.writeKeyword(keyword)
      << frictionContactModelDict_;
  }


// ************************************************************************* //
