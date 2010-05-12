{
    for (label valveI = 0; valveI < nValves(); valveI++)
    {

        // for each valve the following zones have to be created:
        //
        //  - movingFaceZone 
        //  - movingPointZone
        //  - staticFaceZone
        //  - staticPointZone
        
/*
        faceSet movingFaces(*this, valves_[valveI].movingFacesName());
        
        fz.append
        (
            new faceZone
            (
                "movingFacesZoneV" + Foam::name(valveI + 1),
                movingFaces.toc(),
                boolList(movingFaces.size(), false),
                nFaceZones,
                faceZones()
            )
        );
        
        nFaceZones++;
        
        faceSet staticFaces(*this, valves_[valveI].staticFacesName());
        
        fz.append
        (
            new faceZone
            (
                "staticFacesZoneV" + Foam::name(valveI + 1),
                staticFaces.toc(),
                boolList(staticFaces.size(), false),
                nFaceZones,
                faceZones()
            )
        );
        
        nFaceZones++;
*/

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
    }   
}
