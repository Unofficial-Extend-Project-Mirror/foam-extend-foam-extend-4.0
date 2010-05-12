    const labelListList& cp = cellPoints();
    // valves
    forAll(valves_,valveI)
    {
        vector valveVel =
            valves_[valveI].curVelocity()*valves_[valveI].cs().axis();

//         scalar valveDeltaX = valves_[valveI].curVelocity()*
//             valves_[valveI].cs().axis().x()*engTime().deltaT().value() ;
//         scalar valveDeltaY = valves_[valveI].curVelocity()*
//             valves_[valveI].cs().axis().y()*engTime().deltaT().value() ;
        scalar valveDeltaZ = valves_[valveI].curVelocity()*
            valves_[valveI].cs().axis().z()*engTime().deltaT().value() ;

        if(valves_[valveI].curLift() >= valves_[valveI].deformationLift())
        {
            bool isMoving(false);

            if(mag(valves_[valveI].curVelocity()) > 0)
            {
                Info<< "Valve n. " << valveI << " is moving with velocity = "
                    << valves_[valveI].curVelocity() << endl;
                isMoving = true;
            }

            Info<< "Valve displacement for valve " << valveI << " = "
                << valveDeltaZ << endl;

            if
            (
                valves_[valveI].poppetPatchID().active()
             && valves_[valveI].isOpen()
            )
            {
                List<bool> valveTopPoint(newPoints.size(), false);

/*
                label layeringVPtsIndex =
                    pointZones().findZoneID
                    (
                        "movingPointsTopZoneV" + Foam::name(valveI+1)
                    );

                const labelList& layeringVPoints =
                    pointZones()[layeringVPtsIndex];

                forAll(layeringVPoints, i)
                {
                    point& p = newPoints[layeringVPoints[i]];

                    p += valveVel * engTime().deltaT().value();

                    valveTopPoint[layeringVPoints[i]] = true;
                }
*/                

                label layeringVCellsIndex =
                    cellZones().findZoneID
                    (
                        "movingCellsTopZoneV" + Foam::name(valveI+1)
                    );

                const labelList& layeringVCells =
                    cellZones()[layeringVCellsIndex];
                
                labelList layeringVPoints;
                
                boolList count(newPoints.size(), false);
            
                forAll(layeringVCells, cellI)
                {
                    const labelList& curCellPoints = cp[layeringVCells[cellI]];

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

                layeringVPoints.setSize(nCounted);

                // Collect the points
                nCounted = 0;
                forAll (count, pointI)
                {
                    if (count[pointI] == true)
                    {
                        layeringVPoints[nCounted] = pointI;
                        nCounted++;
                    }
                }
                

                forAll(layeringVPoints, i)
                {
                    point& p = newPoints[layeringVPoints[i]];

                    p += valveVel * engTime().deltaT().value();

                    valveTopPoint[layeringVPoints[i]] = true;
                }
                
/*
                label movingTopPtsVIndex =
                    pointZones().findZoneID
                    (
                        "movingPointsV" + Foam::name(valveI+1)
                    );

                const labelList& movingTopPointsV =
                    pointZones()[movingTopPtsVIndex];

                forAll(movingTopPointsV, i)
                {
                    if(!valveTopPoint[movingTopPointsV[i]])
                    {

                        point& p = newPoints[movingTopPointsV[i]];
                        p += valveVel * engTime().deltaT().value();

                    }
                }
*/

                label movingTopCellsVIndex =
                    cellZones().findZoneID
                    (
                        "movingCellsZoneV" + Foam::name(valveI+1)
                    );

                const labelList& movingTopCellsV =
                    cellZones()[movingTopCellsVIndex];

                labelList movingTopPointsV;

                boolList count2(newPoints.size(), false);
            
                forAll(movingTopCellsV, cellI)
                {
                    const labelList& curCellPoints = cp[movingTopCellsV[cellI]];

                    forAll (curCellPoints, i)
                    {
                        count2[curCellPoints[i]] = true;
                    }
                }

                // Count the points
                nCounted = 0;
                forAll (count2, pointI)
                {
                    if (count2[pointI] == true)
                    {
                        nCounted++;
                    }
                }

                movingTopPointsV.setSize(nCounted);

                // Collect the points
                nCounted = 0;
                forAll (count2, pointI)
                {
                    if (count2[pointI] == true)
                    {
                        movingTopPointsV[nCounted] = pointI;
                        nCounted++;
                    }
                }

                forAll(movingTopPointsV, i)
                {
                    if(!valveTopPoint[movingTopPointsV[i]])
                    {

                        point& p = newPoints[movingTopPointsV[i]];
                        p += valveVel * engTime().deltaT().value();

                    }
                }

                
            }
            else
            {
                Info << "Valve " << valveI << " is CLOSED!!!" << endl;
            }

        }
    }
