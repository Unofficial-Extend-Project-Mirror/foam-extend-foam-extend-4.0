    const polyPatch& cylinderHead =
        boundaryMesh()[boundaryMesh().findPatchID("cylinderHead")];

    labelList cylinderHeadFaces(cylinderHead.size());
    forAll(cylinderHeadFaces, i)
    {
        cylinderHeadFaces[i] = cylinderHead.start() + i;
    }

    forAll(valves(), valveI)
    {
        (
            valves_[valveI].detachInCylinderPatchID().active()
         && valves_[valveI].detachInPortPatchID().active()
        )
        {
            Info<< "Adding detach boundary for valve "
                << valveI + 1 << endl;

            const vectorField& areas = Sf().internalField();

            DynamicList<label> df;

            {
                label ndf = 0;
                DynamicList<label> headList;
                label nHead = 0;

                const polyPatch& detachPatch =
                boundaryMesh()[valves()[valveI].poppetPatchID().index()];

                pointField detachPoints(detachPatch.faceCentres());
                detachPoints += valves()[valveI].cs().axis()*valves()[valveI].detachDistance(); 

                forAll(faceCentres(), faceI)
                {
                    if
                    (
                        min(mag(detachPoints - faceCentres()[faceI])) < valves()[valveI].detachTol()
                    )
                    {

                        bool isHead = isACylinderHeadFace(cylinderHeadFaces, faceI);

                        if(isHead)
                        {
                            headList.append(faceI);
                            nHead++;
                        }
                        else
                        {
                            df.append(faceI);
                            ndf++;
                        }
                    }
                }


                Info<< "Found " << ndf << " internal faces and "
                    << nHead << " head faces for the valve n. "
                    << valveI + 1 <<  endl ;

                df.setSize(ndf);
            }

            boolList flip(df.size(), false);

            const vector& valveAxis = valves()[valveI].cs().axis();

            forAll (df, dfI)
            {
                if (isInternalFace(df[dfI]))
                {
                    if ((areas[df[dfI]] &  valveAxis) > 0)
                    {
                        flip[dfI] = true;
                    }
                }
                else
                {
                    FatalErrorIn
                    (
                        "void engineTopoFvMesh::addZonesAndModifiers()"
                    )   << "found boundary face in valve detach definition "
                        << "for valve " << valveI + 1
                        << ".  This is not allowed.  Detach faces: "
                        << df << " nInternalFaces: " << nInternalFaces()
                        << abort(FatalError);
                }
            }

            // Add detach face zone
            fz.append
            (
                new faceZone
                (
                    "detachFaceZoneV" + Foam::name(valveI + 1),
                    df,
                    flip,
                    nFaceZones,
                    faceZones()
                )
            );
            nFaceZones++;
        }
/*


    forAll(valves(), valveI)
    {
        if
        (
            valves_[valveI].detachInCylinderPatchID().active()
         && valves_[valveI].detachInPortPatchID().active()
        )
        {
            Info<< "Adding detach boundary for valve "
                << valveI + 1 << endl;

            const vectorField& areas = Sf().internalField();

            DynamicList<label> df;

            {
                label headFaces = 0;
                label ndf = 0;
                DynamicList<label> headList;
                label nHead = 0;

                const polyPatch& detachPatch =
                boundaryMesh()[valves()[valveI].poppetPatchID().index()];

                forAll(detachPatch.faceCentres(), detachFaceI)
                {
//                    scalar xFacePatch = detachPatch.faceCentres()[detachFaceI].x();
//                    scalar yFacePatch = detachPatch.faceCentres()[detachFaceI].y();
//                    scalar zFacePatch = detachPatch.faceCentres()[detachFaceI].z();

                    point pFacePatch = valves()[valveI].cs().toLocal(detachPatch.faceCentres()[detachFaceI]);

                    scalar xFacePatch = pFacePatch.x();
                    scalar yFacePatch = pFacePatch.y();
                    scalar zFacePatch = pFacePatch.z();
            
                    forAll(faceCentres(),faceI)
                    {
                
                        vector n = faceAreas()[faceI]/mag(faceAreas()[faceI]);


                        scalar dd = n & valves()[valveI].cs().axis();
                        
//                        if (mag(dd) > 0.1)
                        {   
                
                            point pFaceMesh = valves()[valveI].cs().toLocal(faceCentres()[faceI]);

                            scalar xFaceMesh = pFaceMesh.x();
                            scalar yFaceMesh = pFaceMesh.y();
                            scalar zFaceMesh = pFaceMesh.z();
                            
//                            Info <<  "mag(xFaceMesh-xFacePatch) = " << mag(xFaceMesh-xFacePatch) << endl;
                            
                                                 
                            if
                            (
                                mag(xFaceMesh-xFacePatch) < valvePlaneTol_ && 
                                mag(yFaceMesh-yFacePatch) < valvePlaneTol_ && 
                                mag(zFaceMesh-zFacePatch) < 
                                valves()[valveI].detachDistance() + valves()[valveI].detachTol() && 
                                mag(zFaceMesh-zFacePatch) >  
                                valves()[valveI].detachDistance() - valves()[valveI].detachTol()  
//                                && (zFaceMesh-zFacePatch) >  SMALL                          
                            )
                            {
                                bool isHead = isACylinderHeadFace(cylinderHeadFaces, faceI);

                                if(isHead)
                                {
                                    headList.append(faceI);
                                    nHead++;
                                }
                                else
                                {
                                    df.append(faceI);
                                    ndf++;
                                }
                            }
                
                        }
                
                    }
        
                }
                
                Info << "Found " << ndf << " internal faces and " 
                << nHead << " head faces for the valve n. " << valveI + 1 <<  endl ;

                df.setSize(ndf);

            }

            boolList flip(df.size(), false);

            const vector& valveAxis = valves()[valveI].cs().axis();

            forAll (df, dfI)
            {
                if (isInternalFace(df[dfI]))
                {
                    if ((areas[df[dfI]] &  valveAxis) > 0)
                    {
                        flip[dfI] = true;
                    }
                }
                else
                {
                    FatalErrorIn
                    (
                        "void engineTopoFvMesh::addZonesAndModifiers()"
                    )   << "found boundary face in valve detach definition "
                        << "for valve " << valveI + 1
                        << ".  This is not allowed.  Detach faces: "
                        << df << " nInternalFaces: " << nInternalFaces()
                        << abort(FatalError);
                }
            }

            // Add detach face zone
            fz.append
            (
                new faceZone
                (
                    "detachFaceZoneV" + Foam::name(valveI + 1),
                    df,
                    flip,
                    nFaceZones,
                    faceZones()
                )
            );
            nFaceZones++;
        }
        else
        {
            Info << "No valve attach/detach for valve " << valveI + 1 << endl;
        }

*/
    }

    Info << "Added attach/detach faces zones for valves" << endl;

