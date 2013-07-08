//-  the method used to calculate the explicit component of sigma
word divDSigmaLargeStrainExpMethod(mesh.solutionDict().subDict("stressedFoam").lookup("divDSigmaLargeStrainExp"));
Info << "divDSigmaLargeStrainExp calculation method: " << divDSigmaLargeStrainExpMethod << endl;
if(divDSigmaLargeStrainExpMethod != "standard"
   && divDSigmaLargeStrainExpMethod != "surface")
  {
    FatalError << "divDSigmaLargeStrainExp method " << divDSigmaLargeStrainExpMethod << " not found!" << nl
               << "valid methods are:\nstandard\nsurface"
               << exit(FatalError);
  }
