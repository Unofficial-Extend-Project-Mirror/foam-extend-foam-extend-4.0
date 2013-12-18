//- how explicit component of sigma is to be calculated
word sigmaExpMethod(mesh.solutionDict().subDict("solidMechanics").lookup("sigmaExp"));
Info << sigmaExpMethod << " method chosen for calculation of sigmaExp" << endl;
if(sigmaExpMethod != "standard" && sigmaExpMethod != "surface" && sigmaExpMethod != "decompose")
  {
    FatalError << "sigmaExp method " << sigmaExpMethod << " not found!" << nl
        << "valid methods are:\nstandard\nsurface\ndecompose"
        << exit(FatalError);
  }
