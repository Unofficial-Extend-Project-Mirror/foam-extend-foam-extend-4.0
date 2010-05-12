    const polyPatch& cylinderHead =
        boundaryMesh()[boundaryMesh().findPatchID("cylinderHead")];

    labelList cylinderHeadFaces(cylinderHead.size());
    forAll(cylinderHeadFaces, i)
    {
        cylinderHeadFaces[i] = cylinderHead.start() + i;
    }

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
                    scalar xFacePatch = detachPatch.faceCentres()[detachFaceI].x();
                    scalar yFacePatch = detachPatch.faceCentres()[detachFaceI].y();
                    scalar zFacePatch = detachPatch.faceCentres()[detachFaceI].z();

                    forAll(faceCentres(),faceI)
                    {
                        vector n = faceAreas()[faceI]/mag(faceAreas()[faceI]);

                        scalar dd = n & vector(0,0,1);

                        if (mag(dd) > 0.1)
                        {
                            scalar xFaceMesh = faceCentres()[faceI].x();
                            scalar yFaceMesh = faceCentres()[faceI].y();
                            scalar zFaceMesh = faceCentres()[faceI].z();

                            if
                            (
                                mag(xFaceMesh-xFacePatch) < valves()[valveI].detachTol() &&
                                mag(yFaceMesh-yFacePatch) < valves()[valveI].detachTol() &&
                                mag(zFaceMesh-zFacePatch) <
                                valves()[valveI].detachDistance() + valves()[valveI].detachTol() &&
                                mag(zFaceMesh-zFacePatch) >
                                valves()[valveI].detachDistance() - valves()[valveI].detachTol()  &&
                                (zFaceMesh-zFacePatch) >  SMALL
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
                << nHead << " head faces for the valve n. "
                << valveI + 1 <<  endl;

                df.setSize(ndf);
            }

            boolList flip(df.size(), false);

            const vector& pistonAxis = piston().cs().axis();

            forAll (df, dfI)
            {
                if (isInternalFace(df[dfI]))
                {
                    if ((areas[df[dfI]] & pistonAxis) > 0)
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
            fz[nFaceZones] =
                new faceZone
                (
                    "detachFaceZoneV" + Foam::name(valveI + 1),
                    df,
                    flip,
                    nFaceZones,
                    faceZones()
                );
            nFaceZones++;
        }
        else
        {
            Info << "No valve attach/detach for valve " << valveI + 1 << endl;
        }
    }
