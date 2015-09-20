//- write force displacement to file

label leftPatchID = mesh.boundaryMesh().findPatchID("leftClamp");
if(leftPatchID == -1)
  {
    FatalError << "Cannot find patch left for calculating force" << endl;
  }

//- calculate force in x direction on leftClamp patch
scalar leftForce = gSum
(
    vector(1, 0, 0) &
    (mesh.boundary()[leftPatchID].Sf() & sigma.boundaryField()[leftPatchID])
);

//- patchIntegrate utility integrates it this way but this is worng because the sigma tensor should
//- be dotted with the surface normal to give the actual traction/force
//- you cannot just take the component of the sigma tensor
//scalar leftForcePatchIntegrateMethod = gSum
//(
//    mesh.magSf().boundaryField()[leftPatchID]*
//    sigma.boundaryField()[leftPatchID].component(symmTensor::XY)
//);

vector gaugeU1 = vector::zero;
vector gaugeU2 = vector::zero;
if(gaugeFaceID1 != -1)
  {
    gaugeU1 = U.boundaryField()[gaugeFacePatchID1][gaugeFaceID1];
  }
if(gaugeFaceID2 != -1)
  {
    gaugeU2 = U.boundaryField()[gaugeFacePatchID2][gaugeFaceID2];
  }

//- reduce across procs
reduce(gaugeU1, sumOp<vector>());
reduce(gaugeU2, sumOp<vector>());

Pout << "gaugeU1 is " << gaugeU1 << nl
     << "gaugeU2 is " << gaugeU2 << endl;

scalar gaugeDisp = mag(gaugeU1 - gaugeU2);

//- write to file
if(Pstream::master())
  {
    OFstream& forceDispFile = *filePtr;
    forceDispFile << 1000*gaugeDisp << "\t" << -1*leftForce << endl;
  }
