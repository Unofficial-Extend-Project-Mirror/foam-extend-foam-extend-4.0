bool solidInterfaceCorr = rheology.solidInterfaceActive();

solidInterface* solidInterfacePtr(NULL);

if(solidInterfaceCorr)
{
    solidInterfacePtr = &rheology.solInterface();
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
