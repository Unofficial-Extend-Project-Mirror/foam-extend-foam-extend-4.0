/*---------------------------------------------------------------------------*\
correctGlobalFaceZoneMesh.H

When there is a globalFaceZone and the mesh is moved by interpolating U to the
vertices with volPointInterpolation, then there are two problems:
    -some points on the patch with the faceZone are moved incorrectly, I think
     it is because the faceZone has no U and causes an incorrect interpolation,
    -the faceZones points not on the proc cells are not moved at all because
     they have no U.

So the points on the patch with the faceZone need to be fixed and also all the
faceZone points need to be moved and synchronised so each proc has the same
full faceZone mesh.

The mapping of procs faceZone order of points to the master procs faceZone point
order is kept in procToGlobalFZmap, which is calculated at the start of the run
in the createGlobalToLocalFaceZonePointMap.H header file.

Note: DU is used for updated Lagrangian solver instead of U

philipc
\*---------------------------------------------------------------------------*/

//- this is only needed in a parallel runs
if(Pstream::parRun())
  {
    //***** FIX INCORRECT POINT ON PATCHES WITH FACEZONE *****//
    contactPatchPairList& contacts = contact;
    
    forAll(contacts, contactI)
      {
	label masterID = contacts[contactI].masterPatch().index();
	label slaveID = contacts[contactI].slavePatch().index();
	
	primitivePatchInterpolation masterInterpolator
	  (mesh.boundaryMesh()[masterID]);
	primitivePatchInterpolation slaveInterpolator
	  (mesh.boundaryMesh()[slaveID]);

	//- U must be interpolated to the vertices, this ignores the faceZone
	//- points with no U (unlike volPointInterpolation)
	vectorField correctMasterPointU =
	  masterInterpolator.faceToPointInterpolate<vector>
	  (
	   U.boundaryField()[masterID]
	   );
	vectorField correctSlavePointU =
	slaveInterpolator.faceToPointInterpolate<vector>
	  (
	   U.boundaryField()[slaveID]
	   );
	
	vectorField oldMasterPoints =
	  mesh.boundaryMesh()[masterID].localPoints();
	vectorField oldSlavePoints =
	  mesh.boundaryMesh()[slaveID].localPoints();
	
	labelList masterPointLabels =
	  mesh.boundaryMesh()[masterID].meshPoints();
	labelList slavePointLabels =
	  mesh.boundaryMesh()[slaveID].meshPoints();
	
	//- correct the patch newPoints
	forAll(masterPointLabels, pointI)
	  {
	    label pointGlobalLabel = masterPointLabels[pointI];
	    newPoints[pointGlobalLabel] =
	      oldMasterPoints[pointI]
	      +
	      correctMasterPointU[pointI];
	  }
	forAll(slavePointLabels, pointI)
	  {
	    label pointGlobalLabel = slavePointLabels[pointI];
	    newPoints[pointGlobalLabel] =
	      oldSlavePoints[pointI]
	      +
	      correctSlavePointU[pointI];
	  }
      }

    


    //***** NOW FIX AND SYNCHRONISE ALL THE FACEZONE POINTS *****//

    forAll(mesh.faceZones(), faceZoneI)
      {
	//- find the patch corresponding to this faceZone
	//- assuming that the FZ is called <patch_name>FaceZone
	string faceZoneName = mesh.faceZones().names()[faceZoneI];
	//- remove the string FaceZone from the end of the face zone name to get the patch name
	string patchName = faceZoneName.substr(0, (faceZoneName.size()-8));
	label patchID = mesh.boundaryMesh().findPatchID(patchName);
	if(patchID == -1)
	  {
	    FatalError << "Patch " << patchName << " not found corresponding for faceZone"
		       << faceZoneName << exit(FatalError);
	  }
	    
	vectorField globalFZpoints =
	  mesh.faceZones()[faceZoneI]().localPoints();
	
	//- new points for the face zone
	vectorField globalFZnewPoints(globalFZpoints.size(), vector::zero);
	
	//- inter-proc points are shared by multiple procs
	//- pointNumProc is the number of procs which a point lies on
	scalarField pointNumProcs(globalFZpoints.size(), 0.0);
	
	forAll(globalFZnewPoints, globalPointI)
	  {
	    label localPoint = procToGlobalFZmap[faceZoneI][globalPointI];
	    
	    //if(localPoint < mesh.boundaryMesh()[patchID].localPoints().size())
	    if(pointOnLocalProcPatch[faceZoneI][localPoint])
	      {
		label procPoint =
		  mesh.faceZones()[faceZoneI]().meshPoints()[localPoint];
		globalFZnewPoints[globalPointI] =
		  newPoints[procPoint];
		pointNumProcs[globalPointI] = 1;
	      }		
	  }
	
	reduce(globalFZnewPoints, sumOp<vectorField>());
	reduce(pointNumProcs, sumOp<scalarField>());
	
	//- now average the newPoints between all procs
	if(min(pointNumProcs) < 1)
	  {
	    FatalError << "pointNumProc has not been set for all points" << exit(FatalError);
	  }
	globalFZnewPoints /= pointNumProcs;
	
	//- the globalFZnewPoints now contains the correct FZ new points in
	//- a global order, now convert them back into the local proc order
	
	vectorField procFZnewPoints(globalFZpoints.size(), vector::zero);
	
	forAll(globalFZnewPoints, globalPointI)
	  {
	    label localPoint = procToGlobalFZmap[faceZoneI][globalPointI];
	    
	    procFZnewPoints[localPoint] =
	      globalFZnewPoints[globalPointI];
	  }
	
	//- now fix the newPoints points on the globalFaceZones
	labelList procFZmeshPoints =
	  mesh.faceZones()[faceZoneI]().meshPoints();
	
	forAll(procFZmeshPoints, pointI)
	  {
	    label procPoint = procFZmeshPoints[pointI];
	    newPoints[procPoint] =
	      procFZnewPoints[pointI];
	  }
      }
  }
