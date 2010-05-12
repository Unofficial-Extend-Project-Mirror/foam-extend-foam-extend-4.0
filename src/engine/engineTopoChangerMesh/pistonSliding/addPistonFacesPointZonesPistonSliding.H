    // Add the piston zone
    if (piston().patchID().active())
    {

        // Piston position
        
        Info << "Adding face zone for piston layer addition/removal" << endl;
        
        label pistonPatchID = piston().patchID().index();
        
        scalarField zPistonValves(valves().size() + 1,max(boundary()[pistonPatchID].patch().localPoints()).z());
        
        scalarField zPistonValvesOffset = zPistonValves + offSet();
        
        labelListList zonePiston(valves().size()+1);
        boolListList flipZonePiston(valves().size() + 1);
        labelList nZoneFacesPiston(valves().size() + 1, 0);
        boolList foundOneFace(valves().size() + 1, false);        
        scalarList zHigher(valves().size() + 1, GREAT);
        scalarList zLower(valves().size() + 1, GREAT);
        scalarList dh(valves().size() + 1, GREAT);
        scalarList dl(valves().size() + 1, GREAT);
                
        forAll(zonePiston, zoneI)
        {
            zonePiston[zoneI].setSize(faceCentres().size());
            flipZonePiston[zoneI].setSize(faceCentres().size());
            flipZonePiston[zoneI] = false;
        }
        
        forAll (faceCentres(), faceI)
        {
            bool foundInValve = false;
        
            scalar zc = faceCentres()[faceI].z();        
            
            vector n = faceAreas()[faceI]/mag(faceAreas()[faceI]);
            scalar dd = n & vector(0,0,1);
                        
            forAll(valves(), valveI)
            {

                if(inValve(faceCentres()[faceI], valveI))
                {
                    foundInValve = true;
                    
                    scalar zPistV = zPistonValvesOffset[valveI];
                    
                    if (mag(dd) > 0.1)
                    {
                        if (zPistV - zc > 0 && zPistV - zc < dl[valveI])
                        {
                            zLower[valveI] = zc;
                            dl[valveI] = zPistV - zc;
                        }
            
                        if (zc - zPistV > 0 && zc - zPistV < dh[valveI])
                        {
                            zHigher[valveI] = zc;
                            dh[valveI] = zc - zHigher[valveI];
                        }
            
                        if
                        (
                            zc > zPistV - delta()
                            && zc < zPistV + delta()
                        )   
                        {

                            foundOneFace[valveI] = true;

                            if ((faceAreas()[faceI] & vector(0,0,1)) < 0)
                            {
                                flipZonePiston[valveI][nZoneFacesPiston[valveI]] = true;
                            }

                            zonePiston[valveI][nZoneFacesPiston[valveI]] = faceI;
                            nZoneFacesPiston[valveI]++;

                        }
                    }
                }
            }

            if(!foundInValve)
            {
                label valveSize = valves().size();
            
                scalar zPistV = zPistonValvesOffset[valveSize];
                    
                if (mag(dd) > 0.1)
                {
                    if (zPistV - zc > 0 && zPistV - zc < dl[valveSize])
                    {
                        zLower[valveSize] = zc;
                        dl[valveSize] = zPistV - zc;
                    }
            
                    if (zc - zPistV > 0 && zc - zPistV < dh[valveSize])
                    {
                        zHigher[valveSize] = zc;
                        dh[valveSize] = zc - zHigher[valveSize];
                    }
            
                    if
                    (
                        zc > zPistV - delta()
                        && zc < zPistV + delta()
                    )   
                    {

                        foundOneFace[valveSize] = true;

                        if ((faceAreas()[faceI] & vector(0,0,1)) < 0)
                        {
                            flipZonePiston[valveSize][nZoneFacesPiston[valveSize]] = true;
                        }

                        zonePiston[valveSize][nZoneFacesPiston[valveSize]] = faceI;
                        nZoneFacesPiston[valveSize]++;

                    }
                }
            }
        }
        

        // if no cut was found use the layer above

        forAll(valves(), valveI)
        {
        
            if (!foundOneFace[valveI])
            {
                zPistonValvesOffset[valveI] = zHigher[valveI];
            }
        }
        
        if(!foundOneFace[valves().size()])
        {
            zPistonValvesOffset[valves().size()] = zHigher[valves().size()];
        }

        forAll (faceCentres(), faceI)
        {

            bool foundInValve = false;

            scalar zc = faceCentres()[faceI].z();

            vector n = faceAreas()[faceI]/mag(faceAreas()[faceI]);
            scalar dd = n & vector(0,0,1);
                                    
            forAll(valves(), valveI)
            {

                if(!foundOneFace[valveI])
                {
                    
                    scalar zPistV = zPistonValvesOffset[valveI];    
                    
                        
                    if(inValve(faceCentres()[faceI], valveI))
                    {
                        foundInValve = true;
                        
                        if (mag(dd) > 0.1)
                        {   

                            if
                            (
                                zc > zPistV - delta()
                                && zc < zPistV + delta()
                            )
                            {
                                
                                if ((faceAreas()[faceI] & vector(0,0,1)) < 0)
                                {
                                    flipZonePiston[valveI][nZoneFacesPiston[valveI]] = true;
                                }
            
                                zonePiston[valveI][nZoneFacesPiston[valveI]] = faceI;
                                nZoneFacesPiston[valveI]++;
                            }
                        }
                    }

                }
            }
            
            if(!foundInValve && !foundOneFace[valves().size()])
            {
                label valveSize = valves().size();
                scalar zPistV = zPistonValvesOffset[valveSize];    
                
                if (mag(dd) > 0.1)
                {   

                    if
                    (
                        zc > zPistV - delta()
                        && zc < zPistV + delta()
                    )
                    {
                        if ((faceAreas()[faceI] & vector(0,0,1)) < 0)
                        {
                            flipZonePiston[valveSize][nZoneFacesPiston[valveSize]] = true;
                        }
            
                        zonePiston[valveSize][nZoneFacesPiston[valveSize]] = faceI;
                        nZoneFacesPiston[valveSize]++;
                    }
                }
            }
        }

        forAll(valves(), valveI)
        {
            zonePiston[valveI].setSize(nZoneFacesPiston[valveI]);
            flipZonePiston[valveI].setSize(nZoneFacesPiston[valveI]);

            fz.append
            (
                new faceZone
                (
                    "pistonLayerFacesV" + Foam::name(valveI + 1),
                    zonePiston[valveI],
                    flipZonePiston[valveI],
                    nFaceZones,
                    faceZones()
                )    
            );
        
            nFaceZones++;
        }
        
        {
            label valveSize = valves().size();
            zonePiston[valveSize].setSize(nZoneFacesPiston[valveSize]);
            flipZonePiston[valveSize].setSize(nZoneFacesPiston[valveSize]);

            fz.append
            (
                new faceZone
                (
                    "pistonLayerFaces",
                    zonePiston[valveSize],
                    flipZonePiston[valveSize],
                    nFaceZones,
                    faceZones()
                )    
            );
        
            nFaceZones++;
        
        }
        
        forAll(valves(), valveI)        
        {
        
            labelList movingCells(nCells());
                            
            label nMovingCells = 0;
                         
            forAll(cellCentres(),cellI)
            {
                const vector& v = cellCentres()[cellI];
                    
                if
                (
                    inValve(v, valveI)
                    &&
                    v.z() < zPistonValvesOffset[valveI] 
                )
                {
                    movingCells[nMovingCells] = cellI;
                    nMovingCells++;
                }
        
            }
            
            movingCells.setSize(nMovingCells);
            Info << "Number of cells in the moving region poppet valve: " << nMovingCells << endl;
                                            
            cz.append
            (
                new cellZone
                (
                    "movingCellsPistonV"+ Foam::name(valveI + 1),
                    movingCells,
                    nCellZones,
                    cellZones()
                )
            );

            nCellZones++;

            // valve top points (move all with valve displacement)
            DynamicList<label> valvePistonPoints(nPoints() / 10);
            List<bool> valvePistonPoint(nPoints(), false);
            bool foundOne = false;
            
            const cellList& c = cells();
            const faceList& f = allFaces();

            forAll (movingCells, cellI)
            {
                const cell& curCell = c[movingCells[cellI]];

                forAll (curCell, faceI)
                {
                    // Mark all the points as moving
                    const face& curFace = f[curCell[faceI]];

                    forAll (curFace, pointI)
                    {
                        foundOne = true;  
                        valvePistonPoint[curFace[pointI]] = true;            
                    }
        
                }
        
            }
            
            forAll(valvePistonPoint, pointI)
            {
                if(valvePistonPoint[pointI])
                {
                    valvePistonPoints.append(pointI);
                }
            }
                                       
            pz.append
            (
                new pointZone
                (
                    "valvePistonPointsV"+ Foam::name(valveI + 1),
                    valvePistonPoints.shrink(),
                    nPointZones,
                    pointZones()
                )
            );
            
            nPointZones++;

        
        }
        
        
        {
            label valveSize = valves().size();
            labelList movingCells(nCells());
                            
            label nMovingCells = 0;
                         
            forAll(cellCentres(),cellI)
            {
                const vector& v = cellCentres()[cellI];
                
                bool fallInValve = false;
                forAll(valves(), valveI)
                {
                    if(inValve(v, valveI))
                    {
                        fallInValve = true;
                    }
                }
                    
                if
                (
                    !fallInValve
                    &&
                    v.z() < zPistonValvesOffset[valveSize] 
                )
                {
                    movingCells[nMovingCells] = cellI;
                    nMovingCells++;
                }
        
            }
            
            movingCells.setSize(nMovingCells);
            Info << "Number of cells in the moving region poppet valve: " << nMovingCells << endl;
                                            
            cz.append
            (
                new cellZone
                (
                    "movingCellsPiston",
                    movingCells,
                    nCellZones,
                    cellZones()
                )
            );

            nCellZones++;

            // valve top points (move all with valve displacement)
            DynamicList<label> pistonPoints(nPoints() / 10);
            List<bool> pistonPoint(nPoints(), false);
            bool foundOne = false;
            
            const cellList& c = cells();
            const faceList& f = allFaces();

            forAll (movingCells, cellI)
            {
                const cell& curCell = c[movingCells[cellI]];

                forAll (curCell, faceI)
                {
                    // Mark all the points as moving
                    const face& curFace = f[curCell[faceI]];

                    forAll (curFace, pointI)
                    {                            
                        foundOne = true;  
                        pistonPoint[curFace[pointI]] = true;            
                    }
        
                }
        
            }
            
            forAll(pistonPoint, pointI)
            {
                if(pistonPoint[pointI])
                {
                    pistonPoints.append(pointI);
                }
            }
                                        
            pz.append
            (
                new pointZone
                (
                    "pistonPoints",
                    pistonPoints.shrink(),
                    nPointZones,
                    pointZones()
                )
            );
            
            nPointZones++;
        }
        
    }
