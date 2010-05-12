    // valves
    forAll(valves_,valveI)
    {

  	    vector valveVel =
       	    valves_[valveI].curVelocity()*valves_[valveI].cs().axis();
        
        scalar valveDeltaZ = valves_[valveI].curVelocity()*valves_[valveI].cs().axis().z()*engTime().deltaT().value() ;

        if(valves_[valveI].curLift() > valves_[valveI].deformationLift())    
        {        
        

            bool isMoving(false);
        
            if(mag(valves_[valveI].curVelocity()) > 0)
            {
                Info << "Valve n. " << valveI << " is moving with velocity = " << valves_[valveI].curVelocity() << endl;
                isMoving = true;
            }
         
        
            Info << "Valve displacement for valve " << valveI << " = " << valveDeltaZ << endl;

            if(valves_[valveI].poppetPatchID().active() && valves_[valveI].isOpen())
            {

                List<bool> valveTopPoint(newPoints.size(), false);
		
//                label layeringVPtsIndex = pointZones().findZoneID("movingPointsTopZoneV"+ Foam::name(valveI+1));               		
//		        const labelList& layeringVPoints = pointZones()[layeringVPtsIndex];
                
                labelList layeringVPoints;                

                {   
                    label valveLayeringCellsIndex = cellZones().findZoneID("movingCellsTopZoneV" + Foam::name(valveI + 1));
          
                    if (valveLayeringCellsIndex < 0)
                    {      
                        FatalErrorIn("bool verticalValves::update()")
                            << "movingCellsTopZoneV"+ Foam::name(valveI+1)
                            << abort(FatalError);
                    }

                    const labelList& layeringCellsZoneV = cellZones()[valveLayeringCellsIndex];
    
                    const labelListList& cp = cellPoints();

                    boolList count(newPoints.size(), false);

                    forAll (layeringCellsZoneV, cellI)
                    {
                        const labelList& curCellPoints = cp[layeringCellsZoneV[cellI]];

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
              
                }
                
                forAll(layeringVPoints, i)
                {
                    point& p = newPoints[layeringVPoints[i]];
                    p.z() += valveDeltaZ*
                         (valves_[valveI].cs().axis().z()/mag(valves_[valveI].cs().axis()));
                    valveTopPoint[layeringVPoints[i]] = true;
                }
                
                label movingTopPtsVIndex = pointZones().findZoneID("movingPointsV"+ Foam::name(valveI+1));
                		
		        const labelList& movingTopPointsV = pointZones()[movingTopPtsVIndex];
                
                forAll(movingTopPointsV, i)
                {
                    if(!valveTopPoint[movingTopPointsV[i]])
                    {

                        point& p = newPoints[movingTopPointsV[i]];
                        p.z() += valveDeltaZ*
                            (valves_[valveI].cs().axis().z()/mag(valves_[valveI].cs().axis()));
                    }
                }
                		
                deleteDemandDrivenData(movingPointsMaskTopPtr_);

            }
            else
            {
                Info << "Valve " << valveI << " is CLOSED!!!" << endl;
            }
        


       }

       if(valves_[valveI].bottomPatchID().active() )
       {

           {
	
//               label valveBottomPtsIndex = pointZones().findZoneID("movingPointsBottomZoneV"+ Foam::name(valveI+1));          
//	           const labelList& valveBottomPoints = pointZones()[valveBottomPtsIndex];
               
                labelList valveBottomPoints;

                {   
                    label valveBottomCellsIndex = cellZones().findZoneID("movingCellsBottomZoneV" + Foam::name(valveI + 1));
          
                    if (valveBottomCellsIndex < 0)
                    {      
                        FatalErrorIn("bool verticalValves::update()")
                            << "movingCellsBottomZoneV"+ Foam::name(valveI+1)
                            << abort(FatalError);
                    }

                    const labelList& valveMovingCellsBotZoneV = cellZones()[valveBottomCellsIndex];
    
                    const labelListList& cp = cellPoints();

                    boolList count(newPoints.size(), false);

                    forAll (valveMovingCellsBotZoneV, cellI)
                    {
                        const labelList& curCellPoints = cp[valveMovingCellsBotZoneV[cellI]];

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

                    valveBottomPoints.setSize(nCounted);

                    // Collect the points
                    nCounted = 0;
                    forAll (count, pointI)
                    {
                        if (count[pointI] == true)
                        {
                            valveBottomPoints[nCounted] = pointI;
                            nCounted++;
                        }
                    }
              
                }


                forAll(valveBottomPoints, i)
                {
                    point& p = newPoints[valveBottomPoints[i]];
                    p.z() +=    valveDeltaZ*
                           (valves_[valveI].cs().axis().z()/mag(valves_[valveI].cs().axis()));
                }

           }

           deleteDemandDrivenData(movingPointsMaskBottomPtr_);

       }

    }


