{

    if
    (
        piston().patchID().active()
        &&
        piston().bowlInPistonPatchID().active()
        &&
        piston().bowlInCylinderPatchID().active()
    )
    {
        
        Info << "Add sliding interface zones for the piston bowl" << endl;
        
        
        pz.append
        (
            new pointZone
            (
                "cutPointsPistonBowl",
                labelList(0),
                nPointZones,
                pointZones()
            )            
        );
                            
        nPointZones++;

        const polyPatch& bowlInPistonPatch =
            boundaryMesh()
                [piston().bowlInPistonPatchID().index()];

        labelList bowlInPistonLabels(bowlInPistonPatch.size(), bowlInPistonPatch.start());

        forAll (bowlInPistonLabels, i)
        {
            bowlInPistonLabels[i] += i;
        }

        fz.append
        (
            new faceZone
            (
                "bowlInPistonZone",
                bowlInPistonLabels,
                boolList(bowlInPistonLabels.size(), false),
                nFaceZones,
                faceZones()
            )
        );
            
        nFaceZones++;



        const polyPatch& bowlInCylinderPatch =
            boundaryMesh()
                [piston().bowlInCylinderPatchID().index()];

        labelList bowlInCylinderLabels(bowlInCylinderPatch.size(), bowlInCylinderPatch.start());

        forAll (bowlInCylinderLabels, i)
        {
            bowlInCylinderLabels[i] += i;
        }

        fz.append
        (
            new faceZone
            (
                "bowlInCylinderZone",
                bowlInCylinderLabels,
                boolList(bowlInCylinderLabels.size(), false),
                nFaceZones,
                faceZones()
            )
        );
            
        nFaceZones++;



        // Add empty zone for cut faces

        fz.append
        (
            new faceZone
            (
                "cutFaceZonePistonBowl",
                labelList(0),
                boolList(0, false),
                nFaceZones,
                faceZones()
            )    
        );

        nFaceZones++;        
    
    }
    
}
