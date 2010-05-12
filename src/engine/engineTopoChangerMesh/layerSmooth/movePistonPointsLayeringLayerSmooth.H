    Info << "virtualPistonPosition = " << virtualPistonPosition()
    << ", deckHeight = " << deckHeight() << endl;

    // Mesh in three parts:
    // - pistonPoints - move with deltaZ
    // - headPoints - do not move
    
    const pointZoneMesh& pZones = pointZones();
    label headPtsIndex = pZones.findZoneID("headPoints");
    label pistonPtsIndex = pZones.findZoneID("pistonPoints");
    const labelList& pistonPoints = pZones[pistonPtsIndex];
    const labelList& headPoints = pZones[headPtsIndex];


    // Whether point displacement is by scaling
    boolList scaleDisp(nPoints(), true);
    label nScaled = nPoints();
    List<bool> pistonPoint(newPoints.size(), false);
    List<bool> headPoint(newPoints.size(), false);

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
    

    forAll(headPoints, i)
    {
        headPoint[headPoints[i]] = true;
        scaleDisp[headPoints[i]] = false;
        nScaled--;
    }
    
/*
    Info<< "Mesh nPoints:" << nPoints()
        << " inside:" << nScaled
        << " piston:" << pistonPoints.size()
        << " head:" << headPoints.size()
        << endl;
*/

    {

        // Always move piston
        scalar pistonTopZ = -GREAT;
        forAll(pistonPoints, i)
        {
            point& p = newPoints[pistonPoints[i]];
            p.z() += deltaZ;
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
                    if (virtualPistonPosition() > newPoints[i].z())
                    {
                        newPoints[i].z() = pistonTopZ + 0.9*minLayerThickness;
                    }
                }
            }
        }
    }

//    movePoints(newPoints);
//    vpi_.movePoints();
            
    pistonPosition() += deltaZ;
    scalar pistonSpeed = deltaZ/engTime().deltaT().value();

    Info<< "clearance: " << deckHeight() - pistonPosition() << nl
        << "Piston speed = " << pistonSpeed << " m/s" << endl;
