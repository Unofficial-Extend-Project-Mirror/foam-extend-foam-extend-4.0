labelList globalFaceZones;

{
    SLList<label> globalFaceZonesSet;

    const faceZoneMesh& faceZones = stressMesh.faceZones();

    forAll(faceZones, zoneI)
    {
        const faceZone& curFaceZone = faceZones[zoneI];

        forAll(curFaceZone, faceI)
        {
            // if unused face exist
            if (curFaceZone[faceI] >= stressMesh.nFaces())
            {
                globalFaceZonesSet.insert(zoneI);
                break;
            }
        }
    }

    globalFaceZones = labelList(globalFaceZonesSet);
}

labelListList globalToLocalFaceZonePointMap
(
    globalFaceZones.size()
);

forAll(globalFaceZones, zoneI)
{
    label curZoneID = globalFaceZones[zoneI];

    labelList curMap(stressMesh.faceZones()[curZoneID]().nPoints(), -1);

    vectorField fzGlobalPoints = 
        stressMesh.faceZones()[curZoneID]().localPoints();

    //- set all slave points to zero because only the master order is used
    if(!Pstream::master())
    {
        fzGlobalPoints *= 0.0;
    }

    //- pass points to all procs
    reduce(fzGlobalPoints, sumOp<vectorField>());

    //- now every proc has the master's list of FZ points
    //- every proc must now find the mapping from their local FZ points to
    //- the global FZ points

    const vectorField& fzLocalPoints = 
        stressMesh.faceZones()[curZoneID]().localPoints();

    const edgeList& fzLocalEdges =
        stressMesh.faceZones()[curZoneID]().edges();

    const labelListList& fzPointEdges = 
        stressMesh.faceZones()[curZoneID]().pointEdges();

    scalarField minEdgeLength(fzLocalPoints.size(), GREAT);

    forAll(minEdgeLength, pI)
    {
        const labelList& curPointEdges = fzPointEdges[pI];

        forAll(curPointEdges, eI)
        {
            scalar Le = fzLocalEdges[curPointEdges[eI]].mag(fzLocalPoints);
            if (Le < minEdgeLength[pI])
            {
                minEdgeLength[pI] = Le;
            }
        }
    }

    forAll(fzGlobalPoints, globalPointI)
    {
//         scalar minDist = GREAT;

        forAll(fzLocalPoints, procPointI)
        {
            scalar curDist =
                mag
                (
                    fzLocalPoints[procPointI] 
                  - fzGlobalPoints[globalPointI]
                );

//             if (curDist < minDist)
//             {
//                 minDist = curDist;
//             }

            if (curDist < 1e-4*minEdgeLength[procPointI])
            {
                curMap[globalPointI] = procPointI;
                break;
            }
        }

//         if (curMap[globalPointI] == -1)
//         {
//             Pout << "minDist: " << minDist << endl;
//         }
    }

    forAll(curMap, globalPointI)
    {
        if (curMap[globalPointI] == -1)
        {
            FatalErrorIn(args.executable())
                << "local to global face zone point map is not correct" 
                    << abort(FatalError);
        }
    }
    
    globalToLocalFaceZonePointMap[zoneI] = curMap;
}
