//-  the method used to calculate the explicit component of sigma
word divDSigmaLargeStrainExpMethod(mesh.solutionDict().subDict("solidMechanics").lookup("divSigmaLargeStrainExp"));
Info << "divSigmaLargeStrainExp method: " << divDSigmaLargeStrainExpMethod << endl;
if(divDSigmaLargeStrainExpMethod != "standard"
   && divDSigmaLargeStrainExpMethod != "surface")
  {
    FatalError << "divSigmaLargeStrainExp method " << divDSigmaLargeStrainExpMethod << " not found!" << nl
               << "valid methods are:\nstandard\nsurface"
               << exit(FatalError);
  }
