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

            }

        
        }
        
    }
