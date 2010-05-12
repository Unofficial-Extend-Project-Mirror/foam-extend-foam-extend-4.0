//    List<polyMeshModifier*> tm(3*nValves() + 2);
//    DynamicList<polyMeshModifier*> tm;

    label nMods = 0;

    for (label valveI = 0; valveI < nValves(); valveI++)
    {
        
        
        // Add attach-detach for valve
        if
        (
            valves_[valveI].detachInCylinderPatchID().active()
         && valves_[valveI].detachInPortPatchID().active()
        )
        {

            Info << "Adding a attach/detach boundary mesh modifier to the valve " <<
            valveI + 1 << endl;
        
            topoChanger_.setSize(nMods+1);
            topoChanger_.set
            (
                nMods,
                new attachDetach
                (
                    "valveAttachDetach" + Foam::name(valveI + 1),
                    nMods,
                    topoChanger_,
                    "detachFaceZoneV" + Foam::name(valveI + 1),
                    valves_[valveI].detachInCylinderPatchID().name(),
                    valves_[valveI].detachInPortPatchID().name(),
                    scalarField(0),
                    true                // Manual triggering
                )    
            );
            nMods++;
            Info << "Attach detach" << endl;
        }        
                   
    }

    // Add piston layer addition
    if (piston().patchID().active())
    {

/*
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
*/
    }

    Info << "Adding " << nMods << " topology modifiers" << endl;
