bool solidInterfaceCorr = rheology.solidInterfaceActive();

solidInterface* solidInterfacePtr(NULL);

{
    if (solidInterfaceCorr)
    {
        // constitutiveModel is now in charge of solidInterface
        solidInterfacePtr = &rheology.solInterface();
        solidInterfacePtr->modifyProperties(Cf, Kf);

        //- solidInterface needs muf and lambdaf to be used for divSigmaExp
        if (divSigmaExpMethod != "surface")
        {
            FatalError
                << "divSigmaExp must be 'surface' when solidInterface is on"
                << exit(FatalError);
        }

        // check grad scheme
        if (word(mesh.schemesDict().gradSchemes().lookup("grad(U)"))
            != "leastSquaresSolidInterface")
        {
            Warning
                << "The grad(U) gradScheme should be "
                << "leastSquaresSolidInterface for the solidInterface "
                << " procedure to be correct"
                << endl;
        }
    }
}
