    // Add piston layer addition/removal
    if (piston().patchID().active())
    {

        Info << "Adding a layer addition/removal mesh modifier to the piston" << endl;
    
        topoChanger_.setSize(nMods+1);
         
        topoChanger_.set
        (
            nMods,
            new layerAdditionRemoval
            (
                "pistonLayer",
                nMods,
                topoChanger_,
                "pistonLayerFaces",
                piston().minLayer(),
                piston().maxLayer()
            )    
        );
        nMods++;
        Info << "pistonLayer" << endl;
        Info << nMods << endl;            
    }
