{

    const pointZoneMesh& pZones = pointZones();

    label constraintSize = 0;

    forAll(valves_, valveI)
    {
        vector valveVel =
            valves_[valveI].curVelocity()*valves_[valveI].cs().axis();

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

            label movingInternalPtsIndex =
                pZones.findZoneID
                (
                    "movingInternalPointsV"
                  + Foam::name(valveI + 1)
                );

            const labelList& movingInternalPointsV =
                pZones[movingInternalPtsIndex];

            forAll(movingInternalPointsV, mipI)
            {
                constrainedPoints.append(movingInternalPointsV[mipI]);
                constrainedVelocity.append(valveVel);

                constraintSize++;
            }
        }

        {

            label staticPtsIndex =
                pZones.findZoneID("staticPointsV" + Foam::name(valveI + 1));
            const labelList& staticPointsV = pZones[staticPtsIndex];

            forAll(staticPointsV, spI)
            {
                constrainedPoints.append(staticPointsV[spI]);
                constrainedVelocity.append(vector::zero);
                constraintSize++;

            }

        }

        if (piston().patchID().active())
        {
            vector pistonVel =
                piston().cs().axis()*engineTime_.pistonSpeed().value();

            label pistonPtsIndex = pZones.findZoneID("movingPistonPoints");
            const labelList& movingPointsP = pZones[pistonPtsIndex];

            forAll(movingPointsP, mpI)
            {
                constrainedPoints.append(movingPointsP[mpI]);
                constrainedVelocity.append(pistonVel);
//                constrainedVelocity.append(vector::zero);
                constraintSize++;

            }
        }
    }

    Info << "constraintSize = " << constraintSize << endl;
}
