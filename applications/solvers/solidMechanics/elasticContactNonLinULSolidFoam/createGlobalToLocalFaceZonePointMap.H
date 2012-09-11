/*---------------------------------------------------------------------------*\
createGlobalToLocalFaceZonePointMap.H

I could not figure out the order of the points on a globalFaceZone, they is a
different order for each processor where the processors own points come first.

So I decide that the master proc order is the global order and I find the map
from this order to each proc order here by just commparing actual point
coordinates.

This map is then used to correct and synchronise the globalFaceZone points on
all the procs when the mesh is moved, in header file correctGlobalFaceZoneMesh.H.

philipc
\*---------------------------------------------------------------------------*/

//- procToGlobalFZmap is a map from the current proc faceZone point order to the
//- master proc point order
//- these are read if present to allow restarting of contact cases
IOList<labelList> procToGlobalFZmap
(
 IOobject
 (
  "procToGlobalFZmap",
  runTime.timeName(),
  mesh,
  IOobject::READ_IF_PRESENT,
  IOobject::AUTO_WRITE
  ),
 mesh.faceZones().size()
 );

IOList<labelList> pointOnLocalProcPatch
(
 IOobject
 (
  "pointOnLocalProcPatch",
  runTime.timeName(),
  mesh,
  IOobject::READ_IF_PRESENT,
  IOobject::AUTO_WRITE
  ),
 mesh.faceZones().size()
 );

//- if they have been read then don't recalculate it
bool globalFaceZoneMappingSet = false;
if(gMax(procToGlobalFZmap[0]) > 0 && gMax(pointOnLocalProcPatch[0]) > 0)
  {
    Info << "Reading procToGlobalFZmap and pointOnLocalProcPatch allowing restart of contact cases"
	 << endl; 
    globalFaceZoneMappingSet = true;
  }
 else
   {
     Info << "procToGlobalFZmap and pointOnLocalProcPatch will be calculated as it has not been found" << nl
	  << "this message should only appear starting a new analysis" << endl;
   }

//- this is only needed in a parallel runs
if(Pstream::parRun())
  {
    if(!globalFaceZoneMappingSet)
      {
	forAll(mesh.faceZones(), faceZoneI)
	  {
	    vectorField globalFZpoints = mesh.faceZones()[faceZoneI]().localPoints();
	    
	    procToGlobalFZmap[faceZoneI].setSize(globalFZpoints.size(), 0);
	
	    //- set all slave points to zero because only the master order is used
	    if(!Pstream::master())
	      globalFZpoints *= 0.0;
	    
	    //- pass points to all procs
	    reduce(globalFZpoints, sumOp<vectorField>());
	    
	    
	    //- now every proc has the master's list of FZ points
	    //- every proc must now find the mapping from their local FZpoints to
	    //- the globalFZpoints
	    
	    vectorField procFZpoints = mesh.faceZones()[faceZoneI]().localPoints();
	    
	    forAll(globalFZpoints, globalPointI)
	      {
		forAll(procFZpoints, procPointI)
		  {
		    if(procFZpoints[procPointI] == globalFZpoints[globalPointI])
		      {
			procToGlobalFZmap[faceZoneI][globalPointI] = procPointI;
			break;
		      }
		  }
	      }
	    //- procToGlobalFZmap now contains the local FZpoint label for each
	    //- global FZ point label - for each faceZone
	    
	    //- check what points are on the current proc patch
	    pointOnLocalProcPatch[faceZoneI].setSize(globalFZpoints.size(), 0);
	    
	    //- find corresponding patch
	    string faceZoneName = mesh.faceZones().names()[faceZoneI];
	    //- remove the string FaceZone from the end of the face zone name to get the patch name
	    string patchName = faceZoneName.substr(0, (faceZoneName.size()-8));
	    label patchID = mesh.boundaryMesh().findPatchID(patchName);
	    if(patchID == -1)
	      {
		FatalError << "Patch " << patchName << " not found corresponding for faceZone"
			   << faceZoneName << exit(FatalError);
	      }
	    
	    forAll(mesh.faceZones()[faceZoneI]().localPoints(), fzpi)
	      {
		forAll(mesh.boundaryMesh()[patchID].localPoints(), pi)
		  {
		    if(mesh.faceZones()[faceZoneI]().localPoints()[fzpi] == mesh.boundaryMesh()[patchID].localPoints()[pi])
		      {
			pointOnLocalProcPatch[faceZoneI][fzpi] = 1;
			break;
		      }
		  }
	      }
	  }
      } //- end if(!globalFaceZoneMappingSet)
  }

//- write to disk to allow restart of cases
//- because it is not possible to calculate the 
//- mapping after the meshes have moved
if(!globalFaceZoneMappingSet)
  {
    procToGlobalFZmap.write();
    pointOnLocalProcPatch.write();
  }
