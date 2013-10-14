Switch solidInterfaceCorr(false);

solidInterface* solidInterfacePtr(NULL);

{
    const dictionary& stressControl =
        mesh.solutionDict().subDict("solidMechanics");

    solidInterfaceCorr = Switch(stressControl.lookup("solidInterface"));

    if(solidInterfaceCorr)
    {
        Info << "Creating solid interface correction" << endl;
        solidInterfacePtr = new solidInterface(mesh, rheology);
        solidInterfacePtr->modifyProperties(muf, lambdaf);

        //- solidInterface needs muf and lambdaf to be used for divDSigmaExp
        if(divDSigmaExpMethod != "surface" && divDSigmaExpMethod != "decompose")
        {
            FatalError << "divDSigmaExp must be decompose or surface when solidInterface is on"
                << exit(FatalError);
        }
    }
}
