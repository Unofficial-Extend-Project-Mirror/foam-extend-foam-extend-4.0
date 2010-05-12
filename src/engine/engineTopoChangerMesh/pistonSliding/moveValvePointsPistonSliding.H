    // valves
    label constraintSize = 0;

    forAll(valves_,valveI)
    {

          vector valveVel =
               valves_[valveI].curVelocity()*valves_[valveI].cs().axis();

        scalar valveDeltaZ = valves_[valveI].curVelocity()*
            valves_[valveI].cs().axis().z()*engTime().deltaT().value() ;

        if(valves_[valveI].curLift() > valves_[valveI].deformationLift())
        {
            bool isMoving(false);

            if(mag(valves_[valveI].curVelocity()) > 0)
            {
                Info<< "Valve n. " << valveI << " is moving with velocity = "
                    << valves_[valveI].curVelocity() << endl;
                isMoving = true;
            }

            Info<< "Valve displacement for valve " << valveI << " = "
                << valveDeltaZ << endl;

            if
            (
                valves_[valveI].poppetPatchID().active()
             && valves_[valveI].isOpen()
            )
            {
                List<bool> valveTopPoint(newPoints.size(), false);

//                 const scalarField& movingPointsTop =
//                     movingPointsMaskTop(valveI);

                label layeringVPtsIndex =
                    pointZones().findZoneID
                    (
                        "movingPointsTopZoneV"
                      + Foam::name(valveI+1)
                    );

                const labelList& layeringVPoints =
                    pointZones()[layeringVPtsIndex];

                forAll(layeringVPoints, i)
                {
                    point& p = newPoints[layeringVPoints[i]];
                    p.z() += valveDeltaZ*
                         (valves_[valveI].cs().axis().z()/
                         mag(valves_[valveI].cs().axis()));

                    valveTopPoint[layeringVPoints[i]] = true;
                }

                label movingTopPtsVIndex =
                    pointZones().findZoneID
                    (
                        "movingPointsV"
                      + Foam::name(valveI+1)
                    );

                const labelList& movingTopPointsV =
                    pointZones()[movingTopPtsVIndex];

                forAll(movingTopPointsV, i)
                {
                    if(!valveTopPoint[movingTopPointsV[i]])
                    {

                        point& p = newPoints[movingTopPointsV[i]];
                        p.z() += valveDeltaZ*
                            (valves_[valveI].cs().axis().z()/
                            mag(valves_[valveI].cs().axis()));
                    }
                }

                deleteDemandDrivenData(movingPointsMaskTopPtr_);
            }
            else
            {
                Info << "Valve " << valveI << " is CLOSED!!!" << endl;
            }

/*
            {
                label movingPtsIndex =
                    pointZones().findZoneID
                    (
                        "movingPointsV"
                      + Foam::name(valveI + 1)
                    );

                const labelList& movingPointsV = pointZones()[movingPtsIndex];

                forAll(movingPointsV, mpI)
                {
                    constrainedPoints.append(movingPointsV[mpI]);
                    constrainedVelocity.append(valveVel);

                    constraintSize++;
                }
            }
*/
        }

       if (valves_[valveI].bottomPatchID().active())
       {
           const scalarField& movingPointsBottom =
               movingPointsMaskBottom(valveI);

           {
               List<bool> valveBottomPoint(newPoints.size(), false);

               label valveBottomPtsIndex =
                   pointZones().findZoneID
                   (
                       "movingPointsBottomZoneV"
                      + Foam::name(valveI+1)
                   );

               const labelList& valveBottomPoints =
                   pointZones()[valveBottomPtsIndex];

               forAll(valveBottomPoints, i)
               {
                   point& p = newPoints[valveBottomPoints[i]];
                   p.z() +=    valveDeltaZ*
                           (valves_[valveI].cs().axis().z()/
                           mag(valves_[valveI].cs().axis()));
               }

           }

           deleteDemandDrivenData(movingPointsMaskBottomPtr_);
       }
    }
