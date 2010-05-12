
    label nMods = 0;

    for (label valveI = 0; valveI < nValves(); valveI++)
    {
        // Add valve curtain sliding interface
        if
        (
            valves_[valveI].curtainInCylinderPatchID().active()
         && valves_[valveI].curtainInPortPatchID().active()
        )
        {
            
            Info << "Adding a sliding interface mesh modifier to the valve " <<
            valveI + 1 << endl;
            
            topoChanger_.setSize(nMods+1);
            
            topoChanger_.set
            (
                nMods,
                new slidingInterface
                (
                    "valveSlider" + Foam::name(valveI + 1),
                    nMods,
                    topoChanger_,
                    "curtainPortZoneV" + Foam::name(valveI + 1),
                    "curtainCylZoneV" + Foam::name(valveI + 1),
                    "cutPointsV" + Foam::name(valveI + 1),
                    "cutFaceZoneV" + Foam::name(valveI + 1),
                    valves_[valveI].curtainInPortPatchID().name(),
                    valves_[valveI].curtainInCylinderPatchID().name(),
                    slidingInterface::INTEGRAL, // always integral
                    true,  // attach-detach action
                    intersection::VISIBLE
                )    
            );
            nMods++;
            
            Info << "valveSlider" + Foam::name(valveI + 1) << endl;
            Info << nMods << endl;            
        }
        
        // Add valve poppet layer addition
        if (valves_[valveI].poppetPatchID().active())
        {

            Info << "Adding a layer addition/removal mesh modifier to the poppet of the valve " <<
            valveI + 1 << endl;
        
            topoChanger_.setSize(nMods+1);
            
            topoChanger_.set
            (
                nMods,
                new layerAdditionRemoval
                (
                    "valvePoppetLayer" + Foam::name(valveI + 1),
                    nMods,
                    topoChanger_,
                    "poppetZoneV" + Foam::name(valveI + 1),
                    valves_[valveI].minTopLayer(),
                    valves_[valveI].maxTopLayer()
                )    
            );
            nMods++;
            Info << "valvePoppetLayer" + Foam::name(valveI + 1) << endl;
            Info << nMods << endl;            
            
        }

        // Add valve bottom layer addition
        if (valves_[valveI].bottomPatchID().active())
        {
        
            Info << "Adding a layer addition/removal mesh modifier to the bottom of the valve " <<
            valveI + 1 << endl;
        
            topoChanger_.setSize(nMods+1);
            
            topoChanger_.set
            (
                nMods,
                new layerAdditionRemoval
                (
                    "valveBottomLayer" + Foam::name(valveI + 1),
                    nMods,
                    topoChanger_,
                    "bottomZoneV" + Foam::name(valveI + 1),
                    valves_[valveI].minBottomLayer(),
                    valves_[valveI].maxBottomLayer()
                )    
            );
            nMods++;
            Info << "valveBottomLayer" + Foam::name(valveI + 1) << endl;
            Info << nMods << endl;            
            
        }
        
        {

            Info << "Adding a layer addition/removal mesh modifier to the piston-valve " <<
            valveI + 1 << endl;
        
            topoChanger_.setSize(nMods+1);
            
            topoChanger_.set
            (
                nMods,
                new layerAdditionRemoval
                (
                    "valvePistonLayer" + Foam::name(valveI + 1),
                    nMods,
                    topoChanger_,
                    "pistonLayerFacesV" + Foam::name(valveI + 1),
                    piston().minLayer(),
                    piston().maxLayer()
                )    
            );
            nMods++;
            Info << "valvePistonLayer" + Foam::name(valveI + 1) << endl;
            Info << nMods << endl; 
        }   
        
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


    Info << "Adding " << nMods << " topology modifiers" << endl;
