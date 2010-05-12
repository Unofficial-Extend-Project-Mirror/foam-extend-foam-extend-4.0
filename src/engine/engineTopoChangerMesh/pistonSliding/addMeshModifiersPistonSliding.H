//    List<polyMeshModifier*> tm(3*nValves() + 2);
//    DynamicList<polyMeshModifier*> tm;

    label nMods = 0;

    if
    (
        piston().patchID().active()
        &&
        piston().bowlInPistonPatchID().active()
        &&
        piston().bowlInCylinderPatchID().active()
    )
    {
            Info << "Adding a sliding interface mesh modifier to the piston bowl " << endl;
            
            topoChanger_.setSize(nMods+1);
            
            topoChanger_.set
            (
                nMods,
                new slidingInterface
                (
                    "pistonBowlInterface",
                    nMods,
                    topoChanger_,
                    "bowlInCylinderZone",
                    "bowlInPistonZone",
                    "cutPointsPistonBowl",
                    "cutFaceZonePistonBowl",
                    piston().bowlInCylinderPatchID().name(),
                    piston().bowlInPistonPatchID().name(),
                    slidingInterface::INTEGRAL, // always integral
                    true,  // attach-detach action
                    intersection::VISIBLE
                )    
            );
            nMods++;
            
            Info << "pistonBowl" << endl;
            Info << nMods << endl;            
    }

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
                    "layeringFacesTopZoneV" + Foam::name(valveI + 1),
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
                    "layeringFacesBottomZoneV" + Foam::name(valveI + 1),
                    valves_[valveI].minBottomLayer(),
                    valves_[valveI].maxBottomLayer()
                )    
            );
            nMods++;
            Info << "valveBottomLayer" + Foam::name(valveI + 1) << endl;
            Info << nMods << endl;            
            
        }
        
    }	

#   include "addPistonLayerAdditionRemovalMeshModifier.H"
    

    
    Info << "Adding " << nMods << " topology modifiers" << endl;
