//-  the method used to interpolate in cell centre displacements to the vertices to move the mesh
word moveMeshMethod(mesh.solutionDict().subDict("solidMechanics").lookup("moveMeshMethod"));
Info << "Move Mesh method: " << moveMeshMethod << endl;
if(moveMeshMethod != "inverseDistance"
   && moveMeshMethod != "leastSquares")
  {
    FatalError << "moveMeshMethod " << moveMeshMethod << " not found!" << nl
               << "valid methods are:\ninvreseDistance\nleastSquares"
               << exit(FatalError);
  }
