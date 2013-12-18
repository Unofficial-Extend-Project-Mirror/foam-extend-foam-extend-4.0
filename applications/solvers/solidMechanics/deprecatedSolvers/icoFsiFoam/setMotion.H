{
    // Setting mesh motion

    pointVectorField solidPointsDispl = 
        cpi.interpolate(Usolid - Usolid.oldTime());

    vectorField newPoints =
        stressMesh.points() 
      + solidPointsDispl.internalField();

    stressMesh.movePoints(newPoints);


    vectorField fluidPatchPointsDispl =
        interpolatorSolidFluid.pointInterpolate
        (
            solidPointsDispl.boundaryField()[solidPatchID].
            patchInternalField()
        );

    motionUFluidPatch ==
        tppi.pointToPointInterpolate
        (
            fluidPatchPointsDispl/runTime.deltaT().value()
        );

    mesh.update();

#   include "volContinuity.H"

    Info << "Motion magnitude: mean = " 
        << average(mag(Usolid.boundaryField()[solidPatchID]))
        << " max = " 
        << max(mag(Usolid.boundaryField()[solidPatchID])) << endl;
}
