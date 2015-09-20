//- how explicit component of sigma is to be calculated
word divDSigmaNonLinExpMethod(mesh.solutionDict().subDict("solidMechanics").lookup("divSigmaNonLinExp"));
Info << "divSigmaNonLinExp method " << divDSigmaNonLinExpMethod <<  endl;
if(divDSigmaNonLinExpMethod != "standard" && divDSigmaNonLinExpMethod != "surface")
{
    FatalError << "divSigmaNonLinExp method " << divDSigmaNonLinExpMethod << " not found!" << nl
        << "valid methods are:\nstandard\nsurface"
        << exit(FatalError);
}
