{

//    const pointZoneMesh& pZones = pointZones();
    const cellZoneMesh& cZones = cellZones();

    const labelListList& cp = cellPoints();

    label constraintSize = 0;

    forAll(valves_, valveI)
    {
    
        vector valveVel =
            valves_[valveI].curVelocity()*valves_[valveI].cs().axis();

        if(valves_[valveI].curLift() > valves_[valveI].deformationLift())
        {
            valveVel = vector::zero;
        }
    
    
        {
    
//            label movingPtsIndex = pZones.findZoneID("movingPointsV"+Foam::name(valveI + 1));
//            const labelList& movingPointsV = pZones[movingPtsIndex];
/*
            forAll(movingPointsV, mpI)
            {
                constrainedPoints.append(movingPointsV[mpI]);
                constrainedVelocity.append(valveVel);
           
                constraintSize++;
            }
*/
            
            label movingCellsIndex = cZones.findZoneID("movingCellsZoneV"+Foam::name(valveI + 1));
            const labelList& movingCellsV = cZones[movingCellsIndex];
            
            labelList movingPointsV;
            boolList count(newPoints.size(), false);
            
            forAll(movingCellsV, cellI)
            {
                const labelList& curCellPoints = cp[movingCellsV[cellI]];

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

            movingPointsV.setSize(nCounted);

            // Collect the points
            nCounted = 0;
            forAll (count, pointI)
            {
                if (count[pointI] == true)
                {
                    movingPointsV[nCounted] = pointI;
                    nCounted++;
                }
            }

            forAll(movingPointsV, mpI)
            {
                constrainedPoints.append(movingPointsV[mpI]);
                constrainedVelocity.append(valveVel);
           
                constraintSize++;
            }
    
        }

        {
    
//            label staticPtsIndex = pZones.findZoneID("staticPointsV"+Foam::name(valveI + 1));
//            const labelList& staticPointsV = pZones[staticPtsIndex];
/*    
            forAll(staticPointsV, spI)
            {
                constrainedPoints.append(staticPointsV[spI]);
                constrainedVelocity.append(vector::zero);
                constraintSize++;
            }
*/

            label staticCellsIndex = cZones.findZoneID("staticCellsZoneV"+Foam::name(valveI + 1));
            const labelList& staticCellsV = cZones[staticCellsIndex];

            labelList staticPointsV;

            boolList count(newPoints.size(), false);
            
            forAll(staticCellsV, cellI)
            {
                const labelList& curCellPoints = cp[staticCellsV[cellI]];

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

            staticPointsV.setSize(nCounted);

            // Collect the points
            nCounted = 0;
            forAll (count, pointI)
            {
                if (count[pointI] == true)
                {
                    staticPointsV[nCounted] = pointI;
                    nCounted++;
                }
            }
            
            forAll(staticPointsV, spI)
            {
                constrainedPoints.append(staticPointsV[spI]);
                constrainedVelocity.append(vector::zero);
                constraintSize++;
            }

        }

    }
    
    if(piston().patchID().active())
    {
        vector pistonVel =
            piston().cs().axis()*engineTime_.pistonSpeed().value();

/*        
        label pistonPtsIndex = pZones.findZoneID("movingPistonPoints");
        const labelList& movingPointsP = pZones[pistonPtsIndex];
    
        forAll(movingPointsP, mpI)
        {
            constrainedPoints.append(movingPointsP[mpI]);
            constrainedVelocity.append(pistonVel);
            constraintSize++;
        }
*/

        label pistonCellsIndex = cZones.findZoneID("movingPistonCells");
        const labelList& movingCellsP = cZones[pistonCellsIndex];
        
        labelList movingPointsP;
        
        boolList count(newPoints.size(), false);
        
        forAll(movingCellsP, cellI)
        {
            const labelList& curCellPoints = cp[movingCellsP[cellI]];

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

        movingPointsP.setSize(nCounted);

        // Collect the points
        nCounted = 0;
        forAll (count, pointI)
        {
            if (count[pointI] == true)
            {
                movingPointsP[nCounted] = pointI;
                nCounted++;
            }
        }
        
        
        forAll(movingPointsP, mpI)
        {
            constrainedPoints.append(movingPointsP[mpI]);
            constrainedVelocity.append(pistonVel);
            constraintSize++;
        }


    }
    

}
