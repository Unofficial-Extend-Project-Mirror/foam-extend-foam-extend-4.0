//- write force displacement to file

label topPatchID = mesh.boundaryMesh().findPatchID("top");
if(topPatchID == -1)
  {
    FatalError << "Cannot find patch top for calculating force" << endl;
  }

//- calculate force in specified direction on topClamp patch
vector direction(0, 1, 0);

//- for small strain or moving mesh
 scalar topForce = gSum(
 			direction &
 			(mesh.boundary()[topPatchID].Sf() & sigma.boundaryField()[topPatchID])
 			);

//- for large strain total lagrangian
// tensorField F = I + gradU.boundaryField()[topPatchID];
// vectorField totalForce = mesh.Sf().boundaryField()[topPatchID] & (sigma.boundaryField()[topPatchID] & F);

//vector force = sum( totalForce );

//scalar topForce = force[vector::Y];;

//- patchIntegrate utility integrates it this way but this is worng because the sigma tensor should
//- be dotted with the surface normal to give the actual traction/force
//- you cannot just take the component of the sigma tensor
//scalar topForcePatchIntegrateMethod = gSum(
//					    mesh.magSf().boundaryField()[topPatchID]
//					    *sigma.boundaryField()[topPatchID].component(symmTensor::XY)
//					    );

scalar disp = max(U.boundaryField()[topPatchID].component(vector::Y));

//- write to file
if(Pstream::master())
  {
    OFstream& forceDispFile = *filePtr;
    forceDispFile << disp << "\t" << topForce << endl;
  }
