bool solidInterfaceCorr = rheology.solidInterfaceActive();

solidInterface* solidInterfacePtr(NULL);

{
    if (solidInterfaceCorr)
    {
        solidInterfacePtr = &rheology.solInterface();

        //- solidInterface needs muf and lambdaf to be used for divDSigmaExp
        if
        (
            divDSigmaExpMethod != "surface"
         && divDSigmaExpMethod != "decompose"
        )
        {
            FatalErrorIn(args.executable())
                << "divDSigmaExp must be decompose or surface when "
                << "solidInterface is on"
                << exit(FatalError);
        }

        // check grad scheme
        if
        (
            word(mesh.schemesDict().gradSchemes().lookup("grad(U)"))
         != "leastSquaresSolidInterface"
        )
        {
            WarningIn(args.executable())
                << "The grad(U) gradScheme should be "
                << "leastSquaresSolidInterface for the solidInterface "
                << " procedure to be correct"
                << endl;
        }
    }
}
