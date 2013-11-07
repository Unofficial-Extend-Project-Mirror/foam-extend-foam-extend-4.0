bool solidInterfaceCorr = rheology.solidInterfaceActive();

solidInterface* solidInterfacePtr(NULL);

{
    if (solidInterfaceCorr)
    {
        //Info << "Creating solid interface correction" << endl;
        //solidInterfacePtr = new solidInterface(mesh, rheology);
        // constitutiveModel is now in charge of solidInterface
        solidInterfacePtr = &rheology.solInterface();
        solidInterfacePtr->modifyProperties(muf, lambdaf);

        //- solidInterface needs muf and lambdaf to be used for divSigmaExp
        if (divSigmaExpMethod != "surface" && divSigmaExpMethod != "decompose")
        {
            FatalError
                << "divSigmaExp must be decompose or surface when "
                << "solidInterface is on"
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
