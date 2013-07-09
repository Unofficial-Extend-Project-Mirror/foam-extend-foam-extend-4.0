
    // Add the piston points zone
    if (piston().patchID().active())
    {
        IOobject pistonCellsHeader
        (
            piston().pistonCellSetName(),
            time().constant()+"/polyMesh/sets/",
            *this,
            IOobject::MUST_READ
        );

        if(!pistonCellsHeader.headerOk())
        {
            WarningIn
            (
                "accordionEngineMesh::addZonesAndModifiers()"
            )   << "cellSet called " << piston().pistonCellSetName()
                << " does not exist. Continuing mesh motion without piston points" << endl;

            pz.append
            (
                new pointZone
                (
                    "movingPistonPoints",
                    labelList(0),
                    nPointZones,
                    pointZones()
                )
            );

            nPointZones++;

        }
        else
        {
            cellSet movingPistonCells(*this, piston().pistonCellSetName());

            labelList pistonCells = movingPistonCells.toc();

            labelList movingPistonPoints;

            const labelListList& cp = cellPoints();

            boolList count(allPoints().size(), false);

            forAll (pistonCells, cellI)
            {
                const labelList& curCellPoints = cp[pistonCells[cellI]];

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

            movingPistonPoints.setSize(nCounted);

            // Collect the points
            nCounted = 0;
            forAll (count, pointI)
            {
                if (count[pointI] == true)
                {
                    movingPistonPoints[nCounted] = pointI;
                    nCounted++;
                }
            }

            pz.append
            (
                new pointZone
                (
                    "movingPistonPoints",
                    movingPistonPoints,
                    nPointZones,
                    pointZones()
                )
            );

            nPointZones++;

        }

    }

    {

        IOobject headCellsHeader
        (
            headCellSetName_,
            time().constant()+"/polyMesh/sets/",
            *this,
            IOobject::MUST_READ
        );

        if(!headCellsHeader.headerOk())
        {
            WarningIn
            (
                "accordionEngineMesh::addZonesAndModifiers()"
            )   << "cellSet called " << headCellSetName_
                << " does not exist. Continuing mesh motion without head points" << endl;
            pz.append
            (
                new pointZone
                (
                    "headPoints",
                    labelList(0),
                    nPointZones,
                    pointZones()
                )
            );

            nPointZones++;
        }
        else
        {
            cellSet fixedHeadCells(*this, headCellSetName_);

            labelList headCells = fixedHeadCells.toc();

            labelList fixedHeadPoints;

            const labelListList& cp = cellPoints();

            boolList count(allPoints().size(), false);

            forAll (fixedHeadCells, cellI)
            {
                const labelList& curCellPoints = cp[headCells[cellI]];

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

            fixedHeadPoints.setSize(nCounted);

            // Collect the points
            nCounted = 0;
            forAll (count, pointI)
            {
                if (count[pointI] == true)
                {
                    fixedHeadPoints[nCounted] = pointI;
                    nCounted++;
                }
            }

            pz.append
            (
                new pointZone
                (
                    "headPoints",
                    fixedHeadPoints,
                    nPointZones,
                    pointZones()
                )
            );

            nPointZones++;
        }

    }
