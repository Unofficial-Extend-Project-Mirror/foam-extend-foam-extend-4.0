Switch solidInterfaceCorr(false);

solidInterfaceNonLin* solidInterfacePtr(NULL);

{
    const dictionary& stressControl =
        mesh.solutionDict().subDict("solidMechanics");

    solidInterfaceCorr = Switch(stressControl.lookup("solidInterface"));

    if(solidInterfaceCorr)
    {
        Info << "Creating solid interface nonlinear correction" << endl;
        solidInterfacePtr = new solidInterfaceNonLin(mesh, rheology);
        solidInterfacePtr->modifyProperties(muf, lambdaf);

        //- solidInterface needs muf and lambdaf to be used for divDSigmaExp
        if(divDSigmaExpMethod != "surface" && divDSigmaExpMethod != "decompose")
        {
            FatalError << "divDSigmaExp must be decompose or surface when solidInterface is on"
                 << exit(FatalError);
        }
        if(divDSigmaLargeStrainExpMethod != "surface")
        {
            FatalError << "divDSigmaLargeStrainExp must be surface when solidInterface is on"
                << exit(FatalError);
        }
    }
}
