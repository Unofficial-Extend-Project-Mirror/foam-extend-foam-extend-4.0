{
    boolList isConstrained(points().size(), false);

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
                if(!isConstrained[movingPointsV[mpI]])
                {
                    constrainedPoints.append(movingPointsV[mpI]);
                    constrainedVelocity.append(valveVel);

                    constraintSize++;
                    isConstrained[movingPointsV[mpI]] = true;
                }
            }
        }

        {

            label staticPtsIndex =
                pZones.findZoneID("staticPointsV" + Foam::name(valveI + 1));
            const labelList& staticPointsV = pZones[staticPtsIndex];

            forAll(staticPointsV, spI)
            {
                if(!isConstrained[staticPointsV[spI]])
                {
                    constrainedPoints.append(staticPointsV[spI]);
                    constrainedVelocity.append(vector::zero);
                    constraintSize++;
                    isConstrained[staticPointsV[spI]] = true;
                }
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
                if(!isConstrained[movingPointsP[mpI]])
                {
                    constrainedPoints.append(movingPointsP[mpI]);
                    constrainedVelocity.append(pistonVel);
                    constraintSize++;
                    isConstrained[movingPointsP[mpI]] = true;
                }

            }
        }

        {

            label headPtsIndex = pZones.findZoneID("headPoints");
            const labelList& fixedPointsH = pZones[headPtsIndex];

            forAll(fixedPointsH, fpI)
            {
                if(!isConstrained[fixedPointsH[fpI]])
                {
                    constrainedPoints.append(fixedPointsH[fpI]);
                    constrainedVelocity.append(vector::zero);
                    constraintSize++;
                    isConstrained[fixedPointsH[fpI]] = true;
                }

            }
        }
    }

}
