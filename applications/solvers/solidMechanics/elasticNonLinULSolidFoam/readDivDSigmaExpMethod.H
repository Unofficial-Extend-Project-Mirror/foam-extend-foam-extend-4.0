//-  the method used to calculate the explicit component of sigma
word divDSigmaExpMethod(mesh.solutionDict().subDict("solidMechanics").lookup("divSigmaExp"));
Info << "divSigmaExp method: " << divDSigmaExpMethod << endl;
if(divDSigmaExpMethod != "standard"
   && divDSigmaExpMethod != "surface"
   && divDSigmaExpMethod != "decompose"
   && divDSigmaExpMethod != "laplacian")
  {
    FatalError << "divSigmaExp method " << divDSigmaExpMethod << " not found!" << nl
               << "valid methods are:\nstandard\nsurface\ndecompose\nlaplacian"
               << exit(FatalError);
  }
