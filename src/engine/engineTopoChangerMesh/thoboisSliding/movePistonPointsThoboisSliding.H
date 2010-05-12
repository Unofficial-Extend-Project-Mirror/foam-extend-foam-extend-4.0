    {
        
        Info << "moving piston points" << endl;
    
        boolList scaleDisp(nPoints(), true);
        label nScaled = nPoints();

        List<bool> pistonPoint(newPoints.size(), false);

        label pistonPtsIndex = pointZones().findZoneID("pistonPoints");
        const labelList& pistonPoints = pointZones()[pistonPtsIndex];
         
        const scalarField& movingPointsMPiston = movingPointsMaskPiston(); 
         
        forAll(pistonPoints, i)
        {
            label pointI = pistonPoints[i];
            pistonPoint[pointI] = true;
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

            List<bool> pistonPointValve(newPoints.size(), false);

            label pistonValvePtsIndex =
            pointZones().findZoneID("valvePistonPointsV"+Foam::name(valveI + 1));

            const labelList& pistonPointsV = pointZones()[pistonValvePtsIndex];
         
            const scalarField& movingPointsMPistonValves = movingPointsMaskPistonValves(valveI); 

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
