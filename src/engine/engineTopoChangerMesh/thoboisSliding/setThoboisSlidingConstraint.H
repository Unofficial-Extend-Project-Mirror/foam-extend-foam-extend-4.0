{

    const pointZoneMesh& pZones = pointZones();

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
    
//          label movingPtsIndex = pZones.findZoneID("movingPointsV"+Foam::name(valveI + 1));
//          const labelList& movingPointsV = pZones[movingPtsIndex];

            labelList movingPointsV;
          
            {   
                label valveMovingCellsIndex = cellZones().findZoneID("movingCellsZoneV" + Foam::name(valveI + 1));
          
                if (valveMovingCellsIndex < 0)
                {      
                    FatalErrorIn("bool verticalValves::update()")
                        << "movingCellsZoneV"+ Foam::name(valveI+1)
                        << abort(FatalError);
                }

                const labelList& movingCellsZoneV = cellZones()[valveMovingCellsIndex];
    
                const labelListList& cp = cellPoints();

                boolList count(newPoints.size(), false);

                forAll (movingCellsZoneV, cellI)
                {
                    const labelList& curCellPoints = cp[movingCellsZoneV[cellI]];

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
              
            }
    
            forAll(movingPointsV, mpI)
            {
                constrainedPoints.append(movingPointsV[mpI]);
                constrainedVelocity.append(valveVel);
           
                constraintSize++;
            }
    
        }

        {

//        label staticPtsIndex = pZones.findZoneID("staticPointsV"+Foam::name(valveI + 1));
//        const labelList& staticPointsV = pZones[staticPtsIndex];

            labelList staticPointsV;
          
            {   
                label valveStaticCellsIndex = cellZones().findZoneID("staticCellsZoneV" + Foam::name(valveI + 1));
          
                if (valveStaticCellsIndex < 0)
                {      
                    FatalErrorIn("bool verticalValves::update()")
                        << "staticCellsZoneV"+ Foam::name(valveI+1)
                        << abort(FatalError);
                }

                const labelList& staticCellsZoneV = cellZones()[valveStaticCellsIndex];
    
                const labelListList& cp = cellPoints();

                boolList count(newPoints.size(), false);

                forAll (staticCellsZoneV, cellI)
                {
                    const labelList& curCellPoints = cp[staticCellsZoneV[cellI]];

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
              
            }

    
            forAll(staticPointsV, spI)
            {
                constrainedPoints.append(staticPointsV[spI]);
                constrainedVelocity.append(vector::zero);
                constraintSize++;
            }

        }

    }

}
