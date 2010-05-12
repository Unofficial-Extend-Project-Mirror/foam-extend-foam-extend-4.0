    
    // adding face zones for the vertical valves
    
    for (label valveI = 0; valveI < nValves(); valveI++)
    {
        // If both sides of the interface exist, add sliding interface
        // for a valve
        if
        (
            valves_[valveI].curtainInCylinderPatchID().active()
         && valves_[valveI].curtainInPortPatchID().active()
        )
        {
            Info<< "Adding sliding interface zones for curtain of valve "
                << valveI + 1 << endl;

            pz.append
            (
                new pointZone
                (
                    "cutPointsV" + Foam::name(valveI + 1),
                    labelList(0),
                    nPointZones,
                    pointZones()
                )            
            );
                            
            nPointZones++;

            const polyPatch& cylCurtain =
                boundaryMesh()
                    [valves_[valveI].curtainInCylinderPatchID().index()];

            labelList cylCurtainLabels(cylCurtain.size(), cylCurtain.start());

            forAll (cylCurtainLabels, i)
            {
                cylCurtainLabels[i] += i;
            }

            fz.append
            (
                new faceZone
                (
                    "curtainCylZoneV" + Foam::name(valveI + 1),
                    cylCurtainLabels,
                    boolList(cylCurtainLabels.size(), false),
                    nFaceZones,
                    faceZones()
                )
            );
                
            nFaceZones++;

            const polyPatch& portCurtain =
                boundaryMesh()
                    [valves_[valveI].curtainInPortPatchID().index()];

            labelList portCurtainLabels
            (
                portCurtain.size(),
                portCurtain.start()
            );

            forAll (portCurtainLabels, i)
            {
                portCurtainLabels[i] += i;
            }

            fz.append
            (                
                new faceZone
                (
                    "curtainPortZoneV" + Foam::name(valveI + 1),
                    portCurtainLabels,
                    boolList(portCurtainLabels.size(), false),
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
                    "cutFaceZoneV" + Foam::name(valveI + 1),
                    labelList(0),
                    boolList(0, false),
                    nFaceZones,
                    faceZones()
                )    
            );

            nFaceZones++;

        }
        else
        {
            Info << "No valve curtain for valve " << valveI + 1 << endl;
        }
        
     

        // for each valve the following zones have to be created:
        //
        //  - movingFaceZone 
        //  - movingPointZone
        //  - staticFaceZone
        //  - staticPointZone
        

        faceSet layeringFacesTop(*this, valves_[valveI].layeringFacesTopName());
/*
//        if (valves_[valveI].poppetPatchID().active())
//        {
            Info << "Adding face zone for valve top patch layer addition/removal" << endl;

//            label valvePatchID = valves_[valveI].poppetPatchID().index();
                    
//            const polyPatch& valveTopPatch = boundaryMesh()[valvePatchID];
            
//            labelList valveTopLayerFaces(0);
//            boolList flipZone1(0);
            boolList flipZone1(layeringFacesTop.size(), false);
            
            forAll(layeringFacesTop.toc(), faceI)
            {

                    vector n = faceAreas()[faceI]/mag(faceAreas()[faceI]);

                    scalar dd = n & vector(0,0,1);

                    if (mag(dd) > SMALL)
                    {   
                            
                            if ((faceAreas()[faceI] & vector(0,0,1)) > 0)
                            {
                                flipZone1[faceI] = true;
                            }
                            else
                            {
                                flipZone1[faceI] = false;
                            }
                   }
                
                
            }
*/
//        } 
        
	fz.append
        (
            new faceZone
            (
                "layeringFacesTopZoneV" + Foam::name(valveI + 1),
                layeringFacesTop.toc(),
                boolList(layeringFacesTop.size(), true),
                nFaceZones,
                faceZones()
            )
        );
	
        nFaceZones++;

        faceSet layeringFacesBottom(*this, valves_[valveI].layeringFacesBottomName());

        fz.append
        (
            new faceZone
            (
                "layeringFacesBottomZoneV" + Foam::name(valveI + 1),
                layeringFacesBottom.toc(),
                boolList(layeringFacesBottom.size(), true),
                nFaceZones,
                faceZones()
            )
        );
	
        nFaceZones++;

        cellSet movingCellsTop(*this, valves_[valveI].movingCellsTopName());
        
        cz.append
        (
            new cellZone
            (
                "movingCellsTopZoneV" + Foam::name(valveI + 1),
                movingCellsTop.toc(),
                nCellZones,
                cellZones()
            )
        );
        
        nCellZones++;

        cellSet movingCellsBottom(*this, valves_[valveI].movingCellsBottomName());
        
        cz.append
        (
            new cellZone
            (
                "movingCellsBottomZoneV" + Foam::name(valveI + 1),
                movingCellsBottom.toc(),
                nCellZones,
                cellZones()
            )
        );
        
        nCellZones++;

        pointSet movingPointsTop(*this, valves_[valveI].movingPointsTopName());

        pz.append
        (
            new pointZone
            (
                "movingPointsTopZoneV" + Foam::name(valveI + 1),
                movingPointsTop.toc(),
                nPointZones,
                pointZones()
            )
        );
	
        nPointZones++;

        pointSet movingPointsBottom(*this, valves_[valveI].movingPointsBottomName());

        pz.append
        (
            new pointZone
            (
                "movingPointsBottomZoneV" + Foam::name(valveI + 1),
                movingPointsBottom.toc(),
                nPointZones,
                pointZones()
            )
        );
	
        nPointZones++;

        cellSet movingCells(*this, valves_[valveI].movingCellsName());
        
        cz.append
        (
            new cellZone
            (
                "movingCellsZoneV" + Foam::name(valveI + 1),
                movingCells.toc(),
                nCellZones,
                cellZones()
            )
        );
        
        nCellZones++;
        
        cellSet staticCells(*this, valves_[valveI].staticCellsName());
        
        cz.append
        (
            new cellZone
            (
                "staticCellsZoneV" + Foam::name(valveI + 1),
                staticCells.toc(),
                nCellZones,
                cellZones()
            )
        );
        
        nCellZones++;
        

        pointSet movingPoints(*this, valves_[valveI].movingPointsName());

        pz.append
        (
            new pointZone
            (
                "movingPointsV" + Foam::name(valveI + 1),
                movingPoints.toc(),
                nPointZones,
                pointZones()
            )
        );
        
        nPointZones++;        
        
        pointSet staticPoints(*this, valves_[valveI].staticPointsName());

        pz.append
        (
            new pointZone
            (
                "staticPointsV" + Foam::name(valveI + 1),
                staticPoints.toc(),
                nPointZones,
                pointZones()
            )
        );
        
        nPointZones++;     
/*
        pointSet movingInternalPoints(*this, valves_[valveI].movingInternalPointsName());

        pz.append
        (
            new pointZone
            (
                "movingInternalPointsV" + Foam::name(valveI + 1),
                movingInternalPoints.toc(),
                nPointZones,
                pointZones()
            )
        );
        
        nPointZones++;     
*/
    
 }   


