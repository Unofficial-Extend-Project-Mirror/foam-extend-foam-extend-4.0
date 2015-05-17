//- how explicit component of sigma is to be calculated
word divSigmaExpMethod(mesh.solutionDict().subDict("solidMechanics").lookup("divSigmaExp"));
Info << "Selecting divSigmaExp calculation method " << divSigmaExpMethod <<  endl;
if(divSigmaExpMethod != "standard" && divSigmaExpMethod != "surface" && divSigmaExpMethod != "decompose" && divSigmaExpMethod != "laplacian")
  {
    FatalError << "divSigmaExp method " << divSigmaExpMethod << " not found!" << nl
           << "valid methods are:\nstandard\nsurface\ndecompose\nlaplacian"
           << exit(FatalError);
  }
