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

Description
    corrects the contact boundary conditions

\*---------------------------------------------------------------------------*/
#include "contactPatchPair.H"
#include "contactProblem.H"
#include "solidTractionFvPatchVectorField.H"
#include "ListListOps.H"
#include "ggiInterpolation.H"

void Foam::contactPatchPair::correct()
{

  //---------------------PRELIMINARIES---------------------------------//
  const fvMesh& mesh = cp_.mesh();
  const label& masterIndex = masterPatch_.index();
  const label& slaveIndex = slavePatch_.index();
  scalar maxMagSlaveTraction = 0.0;
  contactIterNum_++;


  //--------CALCULATE MASTER AND SLAVE PENETRATIONS----------------------//
  scalarField& globalSlavePointPenetration = globalSlavePointPenetration_;
  //scalarField& globalMasterPointPenetration = globalMasterPointPenetration_;


  //- tell zoneToZone that mesh has moved, so the intersection will be recalculated
  faceZoneMasterToSlaveInterpolator_.movePoints();
  //- calculate intersection distances
  //- this is the slowest part of the contact correction especially when the slavePatch
  //- has many points. parallelisation of this step should be considered.
  globalSlavePointPenetration
    = faceZoneMasterToSlaveInterpolator_.pointDistanceToIntersection();

  //globalMasterPointPenetration
  //= faceZoneSlaveToMasterInterpolator.pointDistanceToIntersection();


  scalarField& slavePointPenetration = slavePointPenetration_;
  //scalarField& masterPointPenetration = masterPointPenetration_;


  forAll(slavePointPenetration, pointI)
    {
      //label pointGlobalLabel = slavePointLabels[pointI];
      slavePointPenetration[pointI] =
	globalSlavePointPenetration
	[
	 pointI //mesh.pointZones()[slavePointZoneID].whichPoint(pointGlobalLabel)
	 ];

      //- when the master surface surrounds the slave (like the pelvis and femur head) then
      //- the slave penetration can sometimes calculate the distance through the femur head
      //- to the pelvis which is wrong so I limit slavePenetration here
      //- i should add a limitPenetration switch here if(limitPenetration)
      if(slavePointPenetration[pointI] < penetrationLimit_)
	{
	  slavePointPenetration[pointI] = 0.0;
	  //globalSlavePointPenetration[pointI] = 0.0;
	}
    }

  //- This is just for visualisation
  // forAll(masterPointPenetration, pointI)
  //   {
  //     masterPointPenetration[pointI] =
  // 	globalMasterPointPenetration
  // 	[
  // 	 pointI
  // 	 ];
  //   }




  //------CALCULATE SLAVE VERTEX FORCES BASED ON PENETRATION-------------//
  //- approximation of penaltyFactor
  //- this should be automatic, these numbers don't really matter, the scaleFactor
  //- scales theses
  scalar bulkModulus = 500e6;
  scalar faceArea = 9e-6; //0.01; //0.0049; //approx
  scalar cellVolume = 2.7e-8; //0.001; //0.00031; //approx
  scalar penaltyFactor = penaltyScaleFactor_*bulkModulus*faceArea*faceArea/cellVolume;
  scalar returnPenaltyFactor = returnScaleFactor_*penaltyFactor;

  //-- slave
  const vectorField& slavePointNormals = mesh.boundaryMesh()[slaveIndex].pointNormals();
  vectorField& totalSlavePointForce = totalSlavePointForce_;

  int numSlaveContactPoints = 0;
  int numSlaveContactPointsReducing = 0;
  int numSlavesUpdated = 0;

  //- so the procs know the global min
  //scalar minSlavePointPenetration = gMin(slavePointPenetration);
  scalar minSlavePointPenetration = gMin(globalSlavePointPenetration);

  {
    //- update old point force
    oldTotalSlavePointForce_ = totalSlavePointForce;

    forAll(totalSlavePointForce, pointI)
      {
	// if a point has penetrated (i.e. if the penetration is negative),
	//add a force to it relative to the penetration
	if(slavePointPenetration[pointI] < -contactGapTol_) //-I had this before < 0.0)
	  {
	    //contactStep = true;
	    numSlaveContactPoints++; // count points in contact
	    numSlavesUpdated++;
	    //- force is linearly dependent on penetration
	    totalSlavePointForce[pointI] +=
	      ( slavePointNormals[pointI] * penaltyFactor * slavePointPenetration[pointI] );
	  }
	//- else if point is within contact tolerance then don't add any more force
	else if(slavePointPenetration[pointI] < 0.0)
	  {
	    numSlaveContactPoints++; // count points in contact
	  }
	// else if penetration is positive and there is a positive
	// pressure (negative traction) still
	// on the point, then slowly reduce the pressure
	else if((totalSlavePointForce[pointI] & slavePointNormals[pointI]) < 0.0)
	  {
	    numSlavesUpdated++;
	    numSlaveContactPointsReducing++;
	    // point forces must be reduced slowly

	    totalSlavePointForce[pointI] +=
	      ( slavePointNormals[pointI] * returnPenaltyFactor * slavePointPenetration[pointI] );

	    // if a tensile force develops
	    if((totalSlavePointForce[pointI] & slavePointNormals[pointI]) > 0.0)
	      {
		totalSlavePointForce[pointI] = vector::zero;
	      }
	  }
      }

    //--------INTERPOLATE SLAVE POINT FORCE TO SLAVE FACES AND APPLY----------//
    //- tell interpolation that mesh has moved
    slaveInterpolator_.movePoints();
    //- local interpolation
    vectorField slaveTraction =
      slaveInterpolator_.pointToFaceInterpolate<vector>
      (
       totalSlavePointForce
       );
    scalarField slavePressure = mag(slaveTraction);

    //- apply slave traction
    {
      solidTractionFvPatchVectorField& slavePatch =
	refCast<solidTractionFvPatchVectorField>(cp_.U().boundaryField()[slaveIndex]);
      slavePatch.traction() = vector::zero;
      slavePatch.pressure() = slavePressure;
      maxMagSlaveTraction = gMax(slavePressure);
    }

    //--------INTERPOLATE SLAVE POINT FORCE TO MASTER FACE TRACTIONS----------//
    //- for a deformable master
    if(!rigidMaster_)
      {
	const label slaveFaceZoneID
	  = mesh.faceZones().findZoneID(slaveFaceZoneName_);
	const label slavePatchStart
	  = mesh.boundaryMesh()[slaveIndex].start();

	scalarField globalSlavePressure
	  (
	   mesh.faceZones()[slaveFaceZoneID].size(),
	   0.0
	   );

	forAll(slavePressure, i)
	  {
	    globalSlavePressure[mesh.faceZones()[slaveFaceZoneID].whichFace(slavePatchStart + i)] =
	      slavePressure[i];
	  }
	//- exchange parallel data
	reduce(globalSlavePressure, maxOp<scalarField>());

	const label masterFaceZoneID = cp_.mesh().faceZones().findZoneID(masterFaceZoneName_);
	scalarField globalMasterPressure(mesh.faceZones()[masterFaceZoneID].size(),0.0);

	if(faceZoneSlaveToMasterInterpolatorPtr_)
	  {
	    zoneToZoneInterpolation& faceZoneSlaveToMasterInterpolator = *faceZoneSlaveToMasterInterpolatorPtr_;
	    faceZoneSlaveToMasterInterpolator.movePoints();
	    //- patchToPatch interpolate tractions - inverse distance weighting
	    globalMasterPressure =
	      faceZoneSlaveToMasterInterpolator.faceInterpolate<scalar>
	      (
	       globalSlavePressure
	       );
	  }
	else if(faceZoneGgiInterpolatorPtr_)
	  {
	    ggiZoneInterpolation& faceZoneGgiInterpolator = *faceZoneGgiInterpolatorPtr_;
	    faceZoneGgiInterpolator.movePoints();

	    //- GGI interpolate tractions
	    globalMasterPressure =
	      faceZoneGgiInterpolator.slaveToMaster
	      (
	       globalSlavePressure
	       );
	  }


	//- exchange parallel data
	reduce(globalMasterPressure, maxOp<scalarField>());

	//Pout << "The max global master trac is " << max(globalMasterPressure) << endl;

	const label masterPatchStart
	  = mesh.boundaryMesh()[masterIndex].start();

	scalarField masterPressure(mesh.boundaryMesh()[masterIndex].size(), 0.0);

	forAll(masterPressure, i)
	  {
	    masterPressure[i] =
	      globalMasterPressure
	      [
	       mesh.faceZones()[masterFaceZoneID].whichFace(masterPatchStart + i)
	       ];
	  }

	//- apply master traction
	{
	  solidTractionFvPatchVectorField& masterPatch =
	    refCast<solidTractionFvPatchVectorField>(cp_.U().boundaryField()[masterIndex]);
	  masterPatch.traction() = vector::zero;
	  masterPatch.pressure() = masterPressure;
	}
      }
    else //- rigid master
      {
	  solidTractionFvPatchVectorField& masterPatch =
	    refCast<solidTractionFvPatchVectorField>(cp_.U().boundaryField()[masterIndex]);
	  masterPatch.traction() = vector::zero;
	  masterPatch.pressure() = 0.0;
      }
  }


  //--------MASTER PROCS WRITES CONTACT INFO FILE----------//
  reduce(numSlaveContactPoints, sumOp<int>());
  reduce(numSlaveContactPointsReducing, sumOp<int>());

  if(Pstream::master())
    {
      OFstream& contactFile = *contactFilePtr_;
//       contactFile << cp_.U().time().value() << "\t\t" << contactStep << "\t\t" << contactIterNum_
// 		  << "\t\t" << penaltyScaleFactor_ << "\t\t" << penaltyFactor << "\t\t" << numSlaveContactPoints
// 		  << "\t\t\t" << minSlavePointPenetration
// 		  << "\t\t" << numSlaveContactPointsReducing << endl;

      int width = 20;
      contactFile << cp_.U().time().value();
      contactFile.width(width);
      contactFile << contactIterNum_;
      contactFile.width(width);
      contactFile << penaltyScaleFactor_;
      contactFile.width(width);
      contactFile << numSlaveContactPoints;
      contactFile.width(width);
      contactFile << minSlavePointPenetration;
      contactFile.width(width);
      contactFile << numSlaveContactPointsReducing;
      contactFile.width(width);
      contactFile << maxMagSlaveTraction << endl;
    }
}
