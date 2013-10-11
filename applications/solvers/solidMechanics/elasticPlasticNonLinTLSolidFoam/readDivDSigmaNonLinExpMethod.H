//- how explicit component of sigma is to be calculated
word divDSigmaNonLinExpMethod(mesh.solutionDict().subDict("stressedFoam").lookup("divDSigmaNonLinExp"));
Info << "Selecting divDSigmaNonLinExp calculation method " << divDSigmaNonLinExpMethod <<  endl;
if(divDSigmaNonLinExpMethod != "standard" && divDSigmaNonLinExpMethod != "surface")
  {
    FatalError << "divDSigmaNonLinExp method " << divDSigmaNonLinExpMethod << " not found!" << nl
	       << "valid methods are:\nstandard\nsurface"
	       << exit(FatalError);
  }
