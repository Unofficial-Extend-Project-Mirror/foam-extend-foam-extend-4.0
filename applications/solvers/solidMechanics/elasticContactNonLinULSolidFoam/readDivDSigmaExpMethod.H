//- how explicit component of sigma is to be calculated
word divDSigmaExpMethod(mesh.solutionDict().subDict("stressedFoam").lookup("divDSigmaExp"));
Info << divDSigmaExpMethod << " method chosen for calculation of DSigmaExp" << endl;
if(divDSigmaExpMethod != "standard" && divDSigmaExpMethod != "surface" && divDSigmaExpMethod != "decompose" && divDSigmaExpMethod != "laplacian")
  {
    FatalError << "divDSigmaExp method " << divDSigmaExpMethod << " not found!" << nl
	       << "valid methods are:\nstandard\nsurface\ndecompose\nlaplacian"
	       << exit(FatalError);
  }
