    {
        // Move fluid mesh
        const vectorField& n = 
            mesh.boundaryMesh()[fluidPatchID].pointNormals();

        primitivePatchInterpolation patchInterpolator
        (
            mesh.boundaryMesh()[fluidPatchID]
        );

        scalarField pointDeltaCoeffs =
            patchInterpolator.faceToPointInterpolate
            (
                mesh.boundary()[fluidPatchID].deltaCoeffs()
            );

        scalar delta = 
            gMax
            (
                mag
                (
                    n
                  & (
                        accumulatedFluidInterfaceDisplacement
                      + fluidPatchPointsDispl
                      - fluidPatchPointsDisplOld
                    )
                )
               *pointDeltaCoeffs
            );

        Info << "Maximal accumulated displacement of interface points: " 
            << delta << endl;

        if(delta < interfaceDeformationLimit)
        {
            // Move only interface points
            pointField newPoints = mesh.allPoints();

            const labelList& meshPoints = 
                mesh.boundaryMesh()[fluidPatchID].meshPoints();

            forAll (fluidPatchPointsDispl, pointI)
            {
                newPoints[meshPoints[pointI]] += 
                    fluidPatchPointsDispl[pointI]
                  - fluidPatchPointsDisplOld[pointI];
            }

            twoDPointCorrector twoDCorrector(mesh);

            twoDCorrector.correctPoints(newPoints);

            mesh.movePoints(newPoints);

            // Accumulate interface points displacement
            accumulatedFluidInterfaceDisplacement += 
                fluidPatchPointsDispl
              - fluidPatchPointsDisplOld;
        }
        else
        {
            // Move whole fluid mesh
            pointField newPoints = mesh.allPoints();

            const labelList& meshPoints = 
                mesh.boundaryMesh()[fluidPatchID].meshPoints();

            forAll (accumulatedFluidInterfaceDisplacement, pointI)
            {
                newPoints[meshPoints[pointI]] -= 
                    accumulatedFluidInterfaceDisplacement[pointI]; 
            }

            twoDPointCorrector twoDCorrector(mesh);

            twoDCorrector.correctPoints(newPoints);

            mesh.movePoints(newPoints);

            accumulatedFluidInterfaceDisplacement += 
                fluidPatchPointsDispl
              - fluidPatchPointsDisplOld;


            if (feMotionSolver)
            {
                tetPointVectorField& motionU =
                    const_cast<tetPointVectorField&>
                    (
                        mesh.objectRegistry::
                        lookupObject<tetPointVectorField>
                        (
                            "motionU"
                        )
                    );
                
                fixedValueTetPolyPatchVectorField& motionUFluidPatch =
                    refCast<fixedValueTetPolyPatchVectorField>
                    (
                        motionU.boundaryField()[fluidPatchID]
                    );

                tetPolyPatchInterpolation tppi
                (
                    refCast<const faceTetPolyPatch>(motionUFluidPatch.patch())
                );

                motionUFluidPatch == 
                    tppi.pointToPointInterpolate
                    (
                        accumulatedFluidInterfaceDisplacement
                       /runTime.deltaT().value()
                    );
            }
            else if (fvMotionSolver)
            {
                pointVectorField& motionU =
                    const_cast<pointVectorField&>
                    (
                        mesh.objectRegistry::
                        lookupObject<pointVectorField>
                        (
                            "pointMotionU"
                        )
                    );
	        
                fixedValuePointPatchVectorField& motionUFluidPatch =
                    refCast<fixedValuePointPatchVectorField>
                    (
                        motionU.boundaryField()[fluidPatchID]
                    );

                motionUFluidPatch == 
                    accumulatedFluidInterfaceDisplacement
                   /runTime.deltaT().value();
            }
            else
            {
                FatalErrorIn(args.executable())
                    << "Problem with mesh motion solver selection"
                        << abort(FatalError);
            }

            mesh.update();

            accumulatedFluidInterfaceDisplacement = vector::zero;
        }
    }
