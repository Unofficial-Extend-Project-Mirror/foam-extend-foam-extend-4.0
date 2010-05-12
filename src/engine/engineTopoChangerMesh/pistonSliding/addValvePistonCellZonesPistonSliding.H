//  Moving mask: moving cells inside the valves....

    {

        for(label i = 1; i <= nValves(); i++)
        {
            
            
            if(valves_[i-1].poppetPatchID().active())
            {
            
                labelList movingCells(nCells());
                                
                label nMovingCells = 0;
        
                label valvePatchID = valves_[i-1].poppetPatchID().index();
                const polyPatch& valveTopPatch = boundaryMesh()[valvePatchID];
                scalar minPopZ = min(boundary()[valvePatchID].patch().localPoints()).z();
                        
                forAll(cellCentres(),cellI)
                {
                    const vector& v = cellCentres()[cellI];
                        
                    if
                    (
                        inValve(v, i-1)
                        &&
                        v.z() > minPopZ 
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
                        "movingCellsTopV"+ Foam::name(i),
                        movingCells,
                        nCellZones,
                        cellZones()
                    )
                );

                nCellZones++;

                // valve top points (move all with valve displacement)
                DynamicList<label> valveTopPoints(nPoints() / 10);
                List<bool> valveTopPoint(nPoints(), false);
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

                            const point& p = points()[curFace[pointI]];
                        
                            forAll(valveTopPatch.localPoints(), j)
                            {
                                scalar xFacePatch = valveTopPatch.localPoints()[j].x();
                                scalar yFacePatch = valveTopPatch.localPoints()[j].y();
                                scalar zFacePatch = valveTopPatch.localPoints()[j].z();
                                
                                if
                                (
                                    mag(p.x() - xFacePatch) < poppetValveTol_
                                    &&
                                    mag(p.y() - yFacePatch) < poppetValveTol_
                                    &&
                                    mag(p.z() - zFacePatch) < valves_[i-1].topLayerOffset() +
                                    valves_[i-1].topLayerTol()
                                )
                                {
                                    foundOne = true;  
                                    valveTopPoint[curFace[pointI]] = true;
                                }
                            }
                
                        }
            
                    }
        
                }
                
                forAll(valveTopPoint, pointI)
                {
                    if(valveTopPoint[pointI])
                    {
                        valveTopPoints.append(pointI);
                    }
                }
                                                
                pz.append
                (
                    new pointZone
                    (
                        "valveTopPointsV"+ Foam::name(i),
                        valveTopPoints.shrink(),
                        nPointZones,
                        pointZones()
                    )
                );
                
                nPointZones++;

            }
            
            if(valves_[i-1].bottomPatchID().active())
            {
            
                labelList movingCells(nCells());
                                
                label nMovingCells = 0;
                
                label valvePatchID = valves_[i-1].bottomPatchID().index();
                const polyPatch& valveBottomPatch = boundaryMesh()[valvePatchID];
                scalar maxPopZ = max(boundary()[valvePatchID].patch().localPoints()).z();


                forAll(cellCentres(),cellI)
                {
                    const vector& v = cellCentres()[cellI];
                        
                    if
                    (
                        inValve(v, i-1)
                        &&
                        v.z() < maxPopZ + SMALL
                    )
                    {
                        movingCells[nMovingCells] = cellI;
                        nMovingCells++;
                    }
            
                }
                
                movingCells.setSize(nMovingCells);
                Info << "Number of cells in the moving region bottom valve: " << nMovingCells << endl;

                cz.append
                (
                    new cellZone
                    (
                        "movingCellsBotV"+ Foam::name(i),
                        movingCells,
                        nCellZones,
                        cellZones()
                    )
                );

                nCellZones++;
                
                
                // valve bottom points (move all with valve displacement)
                DynamicList<label> valveBottomPoints(nPoints() / 10);
                List<bool> valveBottomPoint(nPoints(), false);
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
                            forAll(valveBottomPatch.localPoints(), j)
                            {
                                scalar xFacePatch = valveBottomPatch.localPoints()[j].x();
                                scalar yFacePatch = valveBottomPatch.localPoints()[j].y();
                                scalar zFacePatch = valveBottomPatch.localPoints()[j].z();
                                
                                const point& p = points()[curFace[pointI]];
                                
                                if
                                (
                                    mag(p.x() - xFacePatch) < bottomValveTol_
                                    &&
                                    mag(p.y() - yFacePatch) < bottomValveTol_
                                    &&
                                    mag(p.z() - zFacePatch) < valves_[i-1].bottomLayerOffset() +
                                    valves_[i-1].bottomLayerTol()
                                )
                                {
                                    foundOne = true;  
                                    valveBottomPoint[curFace[pointI]] = true;
                                }
                            }
                
                        }
            
                    }
        
                }

                forAll(valveBottomPoint, pointI)
                {
                    if(valveBottomPoint[pointI])
                    {
                        valveBottomPoints.append(pointI);
                    }
                }
                
                // Add point zone
                                
                pz.append
                (
                    new pointZone
                    (
                        "valveBottomPointsV"+ Foam::name(i),
                        valveBottomPoints.shrink(),
                        nPointZones,
                        pointZones()
                    )
                );

                nPointZones++;

            }

        
        }
        
    }
