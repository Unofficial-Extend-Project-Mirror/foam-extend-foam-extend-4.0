    {
        
        Info << "moving piston points" << endl;
    
        boolList scaleDisp(nPoints(), true);
        label nScaled = nPoints();

        List<bool> pistonPoint(newPoints.size(), false);

        label pistonPtsIndex = pointZones().findZoneID("pistonPoints");
//        const labelList& pistonPoints = pointZones()[pistonPtsIndex];
         
        const scalarField& movingPointsMPiston = movingPointsMaskPiston(); 

        labelList pistonPoints;

        {
            label pistonCellIndex = cellZones().findZoneID("movingCellsPiston");

            if (pistonCellIndex < 0)
            {
                FatalErrorIn("bool verticalValves::update()")
                    << "Cannot find cell zone movingCellsPiston"
                    << abort(FatalError);
            }

            const labelList& pistonCells = cellZones()[pistonCellIndex];

            const labelListList& cp = cellPoints();

            boolList count(newPoints.size(), false);

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

            pistonPoints.setSize(nCounted);

            // Collect the points
            nCounted = 0;
            forAll (count, pointI)
            {
                if (count[pointI] == true)
                {
                    pistonPoints[nCounted] = pointI;
                    nCounted++;
                }
            }

        }

        
           
        {
                        
            // Always move piston
            scalar pistonTopZ = -GREAT;
            forAll(pistonPoints, i)
            {
                point& p = newPoints[pistonPoints[i]];
                p.z() = p.z() + deltaZ*movingPointsMPiston[pistonPoints[i]];
                pistonTopZ = max(pistonTopZ, p.z());
            }

        // NN! fix. only needed for compression
            if (deltaZ > 0.0)
            {
                // check if piston-points have moved beyond the layer above
                forAll(newPoints, i)
                {
                    if (!pistonPoint[i])
                    {
                        bool foundLow = false;
                        if (virtualPistonPosition() > newPoints[i].z())
                        {

                            foundLow = true;
                            
                            newPoints[i].z() 
                            = 
                            (1.0 - movingPointsMPiston[i])*newPoints[i].z() 
                            +
                            movingPointsMPiston[i] *
                            (
                                pistonTopZ 
                            +   
                                (
                                    0.9*minLayerThickness
                                )
                            );
                            
                        }
                    }
                }
            }
        }
        
        deleteDemandDrivenData(movingPointsMaskPistonPtr_);
        
        forAll(valves(), valveI)
        {

            Info << "moving valve points valve n. " << valveI << endl;
        
            boolList scaleDisp(nPoints(), true);
            label nScaled = nPoints();

            const scalarField& movingPointsMPistonValves = movingPointsMaskPistonValves(valveI); 

            List<bool> pistonPointValve(newPoints.size(), false);
            labelList pistonPointsV;

            {
                label pistonCellVIndex = cellZones().findZoneID("movingCellsPistonV"+ Foam::name(valveI + 1));
                
                if (pistonCellVIndex < 0)
                {      
                    FatalErrorIn("bool verticalValves::update()")
                        << "movingCellsPistonV"+ Foam::name(valveI + 1)
                        << abort(FatalError);
                }

                const labelList& pistonCellsV = cellZones()[pistonCellVIndex];

                const labelListList& cp = cellPoints();

                boolList count(newPoints.size(), false);

                forAll (pistonCellsV, cellI)
                {
                    const labelList& curCellPoints = cp[pistonCellsV[cellI]];

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

                pistonPointsV.setSize(nCounted);

                // Collect the points
                nCounted = 0;
                forAll (count, pointI)
                {
                    if (count[pointI] == true)
                    {
                        pistonPointsV[nCounted] = pointI;
                        nCounted++;
                    }
                }
                            
            }


            forAll(pistonPointsV, i)
            {
                label pointI = pistonPointsV[i];
                pistonPointValve[pointI] = true;
                point& p = newPoints[pointI];
            
                if (p.z() < pistonPosition() - 1.0e-6)
                {
                    scaleDisp[pointI] = false;
                    nScaled--;
                }
            }
           
            {
                        
                // Always move piston
                scalar pistonTopZ = -GREAT;
                forAll(pistonPointsV, i)
                {
                    point& p = newPoints[pistonPointsV[i]];
                    
                    p.z() = p.z() + deltaZ*movingPointsMPistonValves[pistonPointsV[i]];
                    pistonTopZ = max(pistonTopZ, p.z());
                }

        // NN! fix. only needed for compression
                if (deltaZ > 0.0)
                {
                    // check if piston-points have moved beyond the layer above
                    forAll(newPoints, i)
                    {
                        if (!pistonPointValve[i])
                        {
                            bool foundLow = false;
                            if (virtualPistonPosition() > newPoints[i].z())
                            {

                                foundLow = true;
                            
                                newPoints[i].z() 
                                = 
                                (1.0 - movingPointsMPistonValves[i])*newPoints[i].z() 
                                +
                                movingPointsMPistonValves[i] *
                                (
                                    pistonTopZ 
                                +      
                                    (
                                        0.9*minLayerThickness
                                    )
                                );
                            
                            }
                        }   
                    }
                }
            }
            
            deleteDemandDrivenData(movingPointsMaskPistonValvesPtr_);
        }
    
    }
