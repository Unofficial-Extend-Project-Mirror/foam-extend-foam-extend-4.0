//- how explicit component of sigma is to be calculated
word divDSigmaExpMethod(mesh.solutionDict().subDict("solidMechanics").lookup("divSigmaExp"));
Info << "Selecting divSigmaExp calculation method " << divDSigmaExpMethod <<  endl;
if(divDSigmaExpMethod != "standard" && divDSigmaExpMethod != "surface" && divDSigmaExpMethod != "decompose" && divDSigmaExpMethod != "laplacian")
{
    FatalError << "divSigmaExp method " << divDSigmaExpMethod << " not found!" << nl
        << "valid methods are:\nstandard\nsurface\ndecompose\nlaplacian"
        << exit(FatalError);
}
