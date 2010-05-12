    
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
        
        // Make a zone for layer addition at the top of the valve
        if (valves_[valveI].poppetPatchID().active())
        {
            Info << "Adding face zone for valve top patch layer addition/removal" << endl;

            label valvePatchID = valves_[valveI].poppetPatchID().index();
                    
            const polyPatch& valveTopPatch = boundaryMesh()[valvePatchID];
            
            labelList valveTopLayerFaces(0);
            boolList flipZone1(0);
            
            forAll(valveTopPatch.faceCentres(), faceI)
            {
                scalar xFacePatch = valveTopPatch.faceCentres()[faceI].x();
                scalar yFacePatch = valveTopPatch.faceCentres()[faceI].y();
                scalar zFacePatch = valveTopPatch.faceCentres()[faceI].z();

                forAll(faceCentres(),faceI)
                {
                
                    vector n = faceAreas()[faceI]/mag(faceAreas()[faceI]);

                    scalar dd = n & vector(0,0,1);

                    if (mag(dd) > SMALL)
                    {   
                        scalar xFaceMesh = faceCentres()[faceI].x();
                        scalar yFaceMesh = faceCentres()[faceI].y();
                        scalar zFaceMesh = faceCentres()[faceI].z();
                
                        if
                        (
                            mag(xFaceMesh-xFacePatch) < poppetValveTol_ && 
                            mag(yFaceMesh-yFacePatch) < poppetValveTol_ && 
                            mag(zFaceMesh-zFacePatch) < valves_[valveI].topLayerOffset() + poppetValveTol_ && 
                            mag(zFaceMesh-zFacePatch) > valves_[valveI].topLayerOffset() - poppetValveTol_  &&
                            (zFaceMesh-zFacePatch) >  SMALL                          
                        )
                        {
                        
                            label sizeTop = valveTopLayerFaces.size();
                            valveTopLayerFaces.setSize(sizeTop + 1);
                            valveTopLayerFaces[sizeTop] = faceI;
                            flipZone1.setSize(sizeTop + 1);
                            if ((faceAreas()[faceI] & vector(0,0,1)) < 0)
                            {
                                flipZone1[sizeTop] = true;
                            }
                            else
                            {
                                flipZone1[sizeTop] = false;
                            }
                                                             
                        }
                    }
                
                }
                
            }

            fz.append
            (
                new faceZone
                (
                    "poppetZoneV" + Foam::name(valveI + 1),
                    valveTopLayerFaces,
                    flipZone1,
                    nFaceZones,
                    faceZones()
                )
            );
        
            nFaceZones++;
                        
        }
        else
        {
            Info << "No poppet layer addition zone for valve "
                << valveI + 1 << endl;
        }
        
        // Make a zone for layer addition at the bottom of the valve
        if (valves_[valveI].bottomPatchID().active())
        {
            Info << "Adding face zone for valve bottom patch layer addition/removal" << endl;

            label valvePatchID = valves_[valveI].bottomPatchID().index();
                    
            const polyPatch& valveBottomPatch = boundaryMesh()[valvePatchID];
            
            labelList valveBottomLayerFaces(0);
            boolList flipZone1(0);
            
            forAll(valveBottomPatch.faceCentres(), faceI)
            {
                scalar xFacePatch = valveBottomPatch.faceCentres()[faceI].x();
                scalar yFacePatch = valveBottomPatch.faceCentres()[faceI].y();
                scalar zFacePatch = valveBottomPatch.faceCentres()[faceI].z();

                forAll(faceCentres(),faceI)
                {
                
                    vector n = faceAreas()[faceI]/mag(faceAreas()[faceI]);

                    scalar dd = n & vector(0,0,1);

                    if (mag(dd) > SMALL)
                    {   
                        scalar xFaceMesh = faceCentres()[faceI].x();
                        scalar yFaceMesh = faceCentres()[faceI].y();
                        scalar zFaceMesh = faceCentres()[faceI].z();
                
                        if
                        (
                            mag(xFaceMesh-xFacePatch) < bottomValveTol_ && 
                            mag(yFaceMesh-yFacePatch) < bottomValveTol_ && 
                            mag(zFaceMesh-zFacePatch) < valves_[valveI].bottomLayerOffset() + bottomValveTol_ && 
                            mag(zFaceMesh-zFacePatch) > valves_[valveI].bottomLayerOffset() - bottomValveTol_  &&
                            (zFacePatch-zFaceMesh) >  SMALL                          
                        )
                        {
                        
                            label sizeTop = valveBottomLayerFaces.size();
                            valveBottomLayerFaces.setSize(sizeTop + 1);
                            valveBottomLayerFaces[sizeTop] = faceI;
                            flipZone1.setSize(sizeTop + 1);
                            if ((faceAreas()[faceI] & vector(0,0,1)) > 0)
                            {
                                flipZone1[sizeTop] = true;
                            }
                            else
                            {
                                flipZone1[sizeTop] = false;
                            }
                                                             
                        }
                    }
                
                }
                
            }

            Info << "valveBottomLayerFaces size = " << valveBottomLayerFaces.size() << endl;    

            fz.append
            (
                new faceZone
                (
                    "bottomZoneV" + Foam::name(valveI + 1),
                    valveBottomLayerFaces,
                    flipZone1,
                    nFaceZones,
                    faceZones()
                )
            );
        
            nFaceZones++;
            
        }
        else
        {
            Info << "No bottom layer addition zone for valve "
                << valveI + 1 << endl;
        }

    }
