//-  the method used to calculate the explicit component of sigma
word divDSigmaExpMethod(mesh.solutionDict().subDict("stressedFoam").lookup("divDSigmaExp"));
Info << "Calculation of divDSigmaExp method: " << divDSigmaExpMethod << endl;
if(divDSigmaExpMethod != "standard"
   && divDSigmaExpMethod != "surface"
   && divDSigmaExpMethod != "decompose"
   && divDSigmaExpMethod != "laplacian")
  {
    FatalError << "divDSigmaExp method " << divDSigmaExpMethod << " not found!" << nl
               << "valid methods are:\nstandard\nsurface\ndecompose\nlaplacian"
               << exit(FatalError);
  }
