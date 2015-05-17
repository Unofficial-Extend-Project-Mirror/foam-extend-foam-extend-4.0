    const dictionary& stressControl =
        mesh.solutionDict().subDict("solidMechanics");

    // Read how explicit component of sigma is to be calculated
    word divDSigmaExpMethod(stressControl.lookup("divSigmaExp"));

    Info<< "Selecting divSigmaExp calculation method: "
        << divDSigmaExpMethod <<  endl;

    if
    (
        divDSigmaExpMethod != "standard"
     && divDSigmaExpMethod != "surface"
     && divDSigmaExpMethod != "decompose"
     && divDSigmaExpMethod != "laplacian"
    )
    {
        FatalErrorIn(args.executable())
            << "divSigmaExp method " << divDSigmaExpMethod << " not found!  "
            << "valid methods are:\nstandard\nsurface\ndecompose\nlaplacian"
            << exit(FatalError);
    }
