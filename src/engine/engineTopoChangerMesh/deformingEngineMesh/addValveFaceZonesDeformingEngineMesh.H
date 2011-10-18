{
    for (label valveI = 0; valveI < nValves(); valveI++)
    {

        // for each valve the following zones have to be created:
        //
        //  - movingFaceZone
        //  - movingPointZone
        //  - staticFaceZone
        //  - staticPointZone

        IOobject movingCellsHeader
        (
            valves_[valveI].movingCellsName(),
            time().constant()+"/polyMesh/sets/",
            *this,
            IOobject::MUST_READ
        );

        if(!movingCellsHeader.headerOk())
        {
            WarningIn
            (
                "deformingEngineMesh::addZonesAndModifiers()"
            )   << "cellSet called " << valves_[valveI].movingCellsName()
                << " does not exist. Continuing mesh motion without rigid motion points for valve " << valves_[valveI].name() << endl;

            pz.append
            (
                new pointZone
                (
                    "movingPointsV" + Foam::name(valveI + 1),
                    labelList(0),
                    nPointZones,
                    pointZones()
                )
            );

            nPointZones++;

        }
        else
        {
            cellSet movingCellsSet(*this, valves_[valveI].movingCellsName());

            labelList movingCells = movingCellsSet.toc();

            labelList movingPoints;

            const labelListList& cp = cellPoints();

            boolList count(allPoints().size(), false);

            forAll (movingCells, cellI)
            {
                const labelList& curCellPoints = cp[movingCells[cellI]];

                forAll (curCellPoints, i)
                {
                    count[curCellPoints[i]] = true;
                }
            }

            // Count the points
            label nCounted = 0;
            forAll (count, pointI)
            {
                if (count[pointI] == true)
                {
                    nCounted++;
                }
            }

            movingPoints.setSize(nCounted);

            // Collect the points
            nCounted = 0;
            forAll (count, pointI)
            {
                if (count[pointI] == true)
                {
                    movingPoints[nCounted] = pointI;
                    nCounted++;
                }
            }

            pz.append
            (
                new pointZone
                (
                    "movingPointsV" + Foam::name(valveI + 1),
                    movingPoints,
                    nPointZones,
                    pointZones()
                )
            );

            nPointZones++;
        }

        IOobject staticCellsHeader
        (
            valves_[valveI].staticCellsName(),
            time().constant()+"/polyMesh/sets/",
            *this,
            IOobject::MUST_READ
        );

        if(!staticCellsHeader.headerOk())
        {
            WarningIn
            (
                "deformingEngineMesh::addZonesAndModifiers()"
            )   << "cellSet called " << valves_[valveI].staticCellsName()
                << " does not exist. Continuing mesh motion without static points for valve " << valves_[valveI].name() << endl;

            pz.append
            (
                new pointZone
                (
                    "staticPointsV" + Foam::name(valveI + 1),
                    labelList(0),
                    nPointZones,
                    pointZones()
                )
            );

            nPointZones++;

        }
        else
        {
            cellSet staticCellSet(*this, valves_[valveI].staticCellsName());

            labelList staticCells = staticCellSet.toc();

            labelList staticPoints;

            const labelListList& cp = cellPoints();

            boolList count(allPoints().size(), false);

            forAll (staticCells, cellI)
            {
                const labelList& curCellPoints = cp[staticCells[cellI]];

                forAll (curCellPoints, i)
                {
                    count[curCellPoints[i]] = true;
                }
            }

            // Count the points
            label nCounted = 0;
            forAll (count, pointI)
            {
                if (count[pointI] == true)
                {
                    nCounted++;
                }
            }

            staticPoints.setSize(nCounted);

            // Collect the points
            nCounted = 0;
            forAll (count, pointI)
            {
                if (count[pointI] == true)
                {
                    staticPoints[nCounted] = pointI;
                    nCounted++;
                }
            }

            pz.append
            (
                new pointZone
                (
                    "staticPointsV" + Foam::name(valveI + 1),
                    staticPoints,
                    nPointZones,
                    pointZones()
                )
            );

            nPointZones++;
        }
    }
}
