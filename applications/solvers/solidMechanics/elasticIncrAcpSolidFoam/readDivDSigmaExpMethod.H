//- how explicit component of sigma is to be calculated
word divDSigmaExpMethod(mesh.solutionDict().subDict("solidMechanics").lookup("divDSigmaExp"));
Info << "Selecting divDSigmaExp calculation method " << divDSigmaExpMethod <<  endl;
if(divDSigmaExpMethod != "standard" && divDSigmaExpMethod != "surface" && divDSigmaExpMethod != "decompose" && divDSigmaExpMethod != "laplacian")
  {
    FatalError << "divDSigmaExp method " << divDSigmaExpMethod << " not found!" << nl
        << "valid methods are:\nstandard\nsurface\ndecompose\nlaplacian"
        << exit(FatalError);
  }
