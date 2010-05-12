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
                
//                if(valves_[valveI].curLift() > valves_[valveI].deformationLift())
//                {
                    const scalarField& movingPointsTop = movingPointsMaskTop(valveI);
                        
                    {

                        List<bool> valveTopPoint(newPoints.size(), false);

                        label valveTopPtsIndex = pointZones().findZoneID("movingPointsTopZoneV"+ Foam::name(valveI+1));
                        const labelList& valveTopPoints = pointZones()[valveTopPtsIndex];

                        forAll(valveTopPoints, i)
                        {
                            label pointI = valveTopPoints[i];
                            valveTopPoint[pointI] = true;
                        }
        
                Info << " valve Delta Z = " << valveDeltaZ << endl;
                        forAll(valveTopPoints, i)
                        {
                            
			    point& p = newPoints[valveTopPoints[i]];
                            p.z() +=    valveDeltaZ*
                                (valves_[valveI].cs().axis().z()/mag(valves_[valveI].cs().axis()))*
                                movingPointsTop[valveTopPoints[i]];
                        }

                    }
            
                    deleteDemandDrivenData(movingPointsMaskTopPtr_);
/*
                }
                else
                {
                    poppetDeformation = true;
                }
*/
            }
            else
            {
                Info << "Valve " << valveI << " is CLOSED!!!" << endl;
            }
        }
    }
