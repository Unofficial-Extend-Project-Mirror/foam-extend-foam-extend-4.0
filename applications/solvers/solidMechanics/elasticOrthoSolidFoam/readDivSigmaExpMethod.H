//- how explicit component of sigma is to be calculated
word divSigmaExpMethod(mesh.solutionDict().subDict("solidMechanics").lookup("divSigmaExp"));
Info << "Calculation of divSigmaExp Method: " << divSigmaExpMethod << endl;
if(divSigmaExpMethod != "standard" && divSigmaExpMethod != "surface" && divSigmaExpMethod != "laplacian")
{
    FatalError << "divSigmaExp method " << divSigmaExpMethod << " not found!" << nl
         << "valid methods are:\nstandard\nsurface\nlaplacian"
         << exit(FatalError);
}
