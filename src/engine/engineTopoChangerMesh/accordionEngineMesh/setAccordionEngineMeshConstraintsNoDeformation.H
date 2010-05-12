{

    const faceZoneMesh& fZones = faceZones();
    const pointZoneMesh& pZones = pointZones();

    label constraintSize = 0;

    forAll(valves_, valveI)
    {
    
        vector valveVel =
            valves_[valveI].curVelocity()*valves_[valveI].cs().axis();

        {
    
            label movingPtsIndex = pZones.findZoneID("movingPointsV"+Foam::name(valveI + 1));
            const labelList& movingPointsV = pZones[movingPtsIndex];
    
            forAll(movingPointsV, mpI)
            {
                constrainedPoints.append(movingPointsV[mpI]);
                constrainedVelocity.append(valveVel);
                constraintSize++;
            }
    
        }

        {

            label staticPtsIndex = pZones.findZoneID("staticPointsV"+Foam::name(valveI + 1));
            const labelList& staticPointsV = pZones[staticPtsIndex];
    
            forAll(staticPointsV, spI)
            {
                constrainedPoints.append(staticPointsV[spI]);
                constrainedVelocity.append(vector::zero);
                constraintSize++;
            }

        }

        {
    
            label movingFacesIndex = fZones.findZoneID("movingFacesZoneV"+Foam::name(valveI + 1));
            const labelList& movingFacesV = pZones[movingFacesIndex];
    
            forAll(movingFacesV, mfI)
            {
                constrainedPoints.append(movingFacesV[mfI]+points().size());
                constrainedVelocity.append(valveVel);
                constraintSize++;
            }

        }

        {
    
            label staticFacesIndex = fZones.findZoneID("staticFacesZoneV"+Foam::name(valveI + 1));
            const labelList& staticFacesV = pZones[staticFacesIndex];
    
            forAll(staticFacesV, sfI)
            {
                constrainedPoints.append(staticFacesV[sfI]+points().size());
                constrainedVelocity.append(vector::zero);
                constraintSize++;
            }

        }

    }

}
