//- net forces on boundary patches

vectorField netForces(mesh.boundary().size(), vector::zero);

Info << nl;
forAll(netForces, patchI)
{
  netForces[patchI] = gSum(mesh.Sf().boundaryField()[patchI] & sigma.boundaryField()[patchI]);

  Info << "patch\t" << mesh.boundary()[patchI].name() << "\t\tnet force is\t"
       << netForces[patchI] << " N" << endl;
}
