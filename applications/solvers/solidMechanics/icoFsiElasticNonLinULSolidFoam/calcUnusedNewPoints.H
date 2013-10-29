forAll(globalFaceZones, zoneI)
{
    const label curZoneID = globalFaceZones[zoneI];

    const labelList& curMap = 
        globalToLocalFaceZonePointMap[zoneI];

    const labelList& curZoneMeshPoints =
        stressMesh.faceZones()[curZoneID]().meshPoints();

    vectorField curGlobalZonePointDispl
    (
        curZoneMeshPoints.size(), 
        vector::zero
    );

    //-Inter-proc points are shared by multiple procs
    // pointNumProc is the number of procs which a point lies on
    scalarField pointNumProcs(curZoneMeshPoints.size(), 0);

    forAll(curGlobalZonePointDispl, globalPointI)
    {
        label localPoint = curMap[globalPointI];
	    
        if(curZoneMeshPoints[localPoint] < stressMesh.nPoints())
        {
            label procPoint = curZoneMeshPoints[localPoint];
            
            curGlobalZonePointDispl[globalPointI] = pointDUI[procPoint];
            
            pointNumProcs[globalPointI] = 1;
        }
    }

    if (Pstream::parRun())
    {
        reduce(curGlobalZonePointDispl, sumOp<vectorField>());
        reduce(pointNumProcs, sumOp<scalarField>());

        //- now average the displacement between all procs
        curGlobalZonePointDispl /= pointNumProcs;
    }

    //- The curZonePointsDisplGlobal now contains the correct face zone 
    //  displacement in a global master processor order, now convert them 
    //  back into the local proc order

    vectorField curZonePointDispl
    (
        curZoneMeshPoints.size(), 
        vector::zero
    );

    forAll(curGlobalZonePointDispl, globalPointI)
    {
        label localPoint = curMap[globalPointI];

        curZonePointDispl[localPoint] = 
            curGlobalZonePointDispl[globalPointI];
    }

    forAll(curZonePointDispl, pointI)
    {
        // unused points
        if (curZoneMeshPoints[pointI] >= stressMesh.nPoints())
        {
            newPoints[curZoneMeshPoints[pointI]] +=
                curZonePointDispl[pointI];
        }
    }
}
