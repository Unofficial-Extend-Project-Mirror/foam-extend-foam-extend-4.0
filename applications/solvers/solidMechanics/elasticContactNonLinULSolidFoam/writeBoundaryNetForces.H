// * * * * * * * * * * * * * * * * NET FORCES * * * * * * * * * * * * * * *  //

vectorField netForces(mesh.boundary().size(), vector::zero);

Info << nl;
forAll(netForces, patchI)
{
  netForces[patchI] = gSum(mesh.Sf().boundaryField()[patchI] & sigma.boundaryField()[patchI]);
  
  Info << "patch " << mesh.boundary()[patchI].name() << " net force is "
       << netForces[patchI] << " N" << endl;
}   

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
