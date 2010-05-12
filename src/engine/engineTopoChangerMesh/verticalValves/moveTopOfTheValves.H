    {
        forAll(valves(), valveI)
        {
            scalar valveDeltaZ = valves_[valveI].curVelocity()*engTime().deltaT().value() ;
            bool isMoving(false);
        
            if(mag(valves_[valveI].curVelocity()) > 0)
            {
                Info << "Valve n. " << valveI << " is moving" << endl;
                isMoving = true;
            }
            
            if(valves_[valveI].poppetPatchID().active() && valves_[valveI].isOpen())
            {
                
                if(valves_[valveI].curLift() > valves_[valveI].deformationLift())
                {
                    const scalarField& movingPointsTop = movingPointsMaskTop(valveI);
                        
                    {

                        List<bool> valveTopPoint(newPoints.size(), false);

                        labelList valveTopPoints;
                
                        {   
                            label valveTopCellsIndex = cellZones().findZoneID("movingCellsTopV"+ Foam::name(valveI+1));
                
                            if (valveTopCellsIndex < 0)
                            {      
                                FatalErrorIn("bool verticalValves::update()")
                                    << "movingCellsTopV"+ Foam::name(valveI+1)
                                    << abort(FatalError);
                            }

                            const labelList& movingCellsTopV = cellZones()[valveTopCellsIndex];
    
                            const labelListList& cp = cellPoints();

                            boolList count(newPoints.size(), false);

                            forAll (movingCellsTopV, cellI)
                            {
                                const labelList& curCellPoints = cp[movingCellsTopV[cellI]];

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

                            valveTopPoints.setSize(nCounted);

                            // Collect the points
                            nCounted = 0;
                            forAll (count, pointI)
                            {
                                if (count[pointI] == true)
                                {
                                    valveTopPoints[nCounted] = pointI;
                                    nCounted++;
                                }
                            }
                            
                        }


//                        label valveTopPtsIndex = pointZones().findZoneID("valveTopPointsV"+ Foam::name(valveI+1));
//                        const labelList& valveTopPoints = pointZones()[valveTopPtsIndex];

                        forAll(valveTopPoints, i)
                        {
                            label pointI = valveTopPoints[i];
                            valveTopPoint[pointI] = true;
                        }
        
                
                        forAll(valveTopPoints, i)
                        {
                            point& p = newPoints[valveTopPoints[i]];
                            p.z() +=    valveDeltaZ*
                                (valves_[valveI].cs().axis().z()/mag(valves_[valveI].cs().axis()))*
                                movingPointsTop[valveTopPoints[i]];
                        }

                    }
            
                    deleteDemandDrivenData(movingPointsMaskTopPtr_);
                }
                else
                {
                    poppetDeformation = true;
                }
            }
            else
            {
                Info << "Valve " << valveI << " is CLOSED!!!" << endl;
            }
        }
    }
