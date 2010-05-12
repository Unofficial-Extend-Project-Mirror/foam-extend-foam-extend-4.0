{

    const faceZoneMesh& fZones = faceZones();
    const pointZoneMesh& pZones = pointZones();
    const cellZoneMesh& cZones = cellZones();

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
            label movingPtsIndex =
                pZones.findZoneID("movingPointsV"+Foam::name(valveI + 1));

            const labelList& movingPointsV = pZones[movingPtsIndex];

            forAll(movingPointsV, mpI)
            {
                constrainedPoints.append(movingPointsV[mpI]);
                constrainedVelocity.append(valveVel);

                constraintSize++;
            }
        }

        {
            label staticPtsIndex =
                pZones.findZoneID("staticPointsV"+Foam::name(valveI + 1));

            const labelList& staticPointsV = pZones[staticPtsIndex];

            forAll(staticPointsV, spI)
            {
                constrainedPoints.append(staticPointsV[spI]);
                constrainedVelocity.append(vector::zero);
                constraintSize++;
            }
        }

        {

/*
            label movingBottomPtsIndex =
                pZones.findZoneID
                (
                    "movingPointsBottomZoneV"
                   + Foam::name(valveI + 1)
                );

            const labelList& movingBottomPointsV =
                pZones[movingBottomPtsIndex];

            forAll(movingBottomPointsV, mbpI)
            {
                constrainedPoints.append(movingBottomPointsV[mbpI]);
                constrainedVelocity.append(valveVel);
                constraintSize++;
            }
*/
        }
    }
}
