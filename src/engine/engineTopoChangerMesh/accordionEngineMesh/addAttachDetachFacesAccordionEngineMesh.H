
    for (label valveI = 0; valveI < nValves(); valveI++)
    {
        
        label nFlipDetach = 0;
        
        faceSet detachFaces(*this, valves_[valveI].detachFacesName());
        
        boolList flipDetach(detachFaces.size(), false);
                   
        const polyPatch& bottomPatch =
        boundaryMesh()[valves()[valveI].bottomPatchID().index()];
        
        vector bottomPatchCenter = sum(mag(bottomPatch.faceAreas()) *
        bottomPatch.faceCentres())/sum(mag(bottomPatch.faceAreas()));
                
        forAll(detachFaces.toc(), i)
        {
            vector radius = faceCentres()[detachFaces.toc()[i]] - bottomPatchCenter;
            
            radius/= mag(radius);
            
            if(((faceAreas()[detachFaces.toc()[i]]) & (radius) / mag(faceAreas()[detachFaces.toc()[i]])) < 0)
            {
//                flipDetach[i] = true;
//                nFlipDetach++;
            }
               
        }
        
        Info << "valve n. " << valveI << " found " << nFlipDetach << " flipped faces for detach" << endl;
                
        fz.append
        (
            new faceZone
            (
                "detachFaceZoneV" + Foam::name(valveI + 1),
                detachFaces.toc(),
                flipDetach,
                nFaceZones,
                faceZones()
            )
        );
        
        nFaceZones++;     
    }   

    Info << "Added attach/detach faces zones for valves" << endl;
    
