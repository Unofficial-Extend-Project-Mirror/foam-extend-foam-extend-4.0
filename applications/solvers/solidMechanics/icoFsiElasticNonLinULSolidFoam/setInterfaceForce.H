{
    Info << "Setting traction on solid patch" << endl;

//     vectorField fluidPatchTraction =
//        -rhoFluid.value()*nu.value()
//        *U.boundaryField()[fluidPatchID].snGrad()
//       + rhoFluid.value()*p.boundaryField()[fluidPatchID]
//        *mesh.boundary()[fluidPatchID].nf();

    vectorField fluidPatchTraction =
       -rhoFluid.value()*nu.value()
       *U.boundaryField()[fluidPatchID].snGrad();

    scalarField fluidPatchPressure =
        rhoFluid.value()*p.boundaryField()[fluidPatchID];

    vectorField fluidZoneTraction
        (
            mesh.faceZones()[fluidZoneID].size(), 
            vector::zero
        );

    const label fluidPatchStart = 
        mesh.boundaryMesh()[fluidPatchID].start();

    forAll(fluidPatchTraction, i)
    {
        fluidZoneTraction
        [
            mesh.faceZones()[fluidZoneID].whichFace(fluidPatchStart + i)
        ] = 
            fluidPatchTraction[i];
    }

    // Parallel data exchange: collect pressure field on all processors
    reduce(fluidZoneTraction, sumOp<vectorField>());


    scalarField fluidZonePressure
        (
            mesh.faceZones()[fluidZoneID].size(), 
            0.0
        );

    forAll(fluidPatchPressure, i)
    {
        fluidZonePressure
        [
            mesh.faceZones()[fluidZoneID].whichFace(fluidPatchStart + i)
        ] = 
            fluidPatchPressure[i];
    }

    // Parallel data exchange: collect pressure field on all processors
    reduce(fluidZonePressure, sumOp<scalarField>());

    vectorField solidZoneTraction =
        interpolatorFluidSolid.faceInterpolate
        (
            fluidZoneTraction
        );

    scalarField solidZonePressure =
        interpolatorFluidSolid.faceInterpolate
        (
            fluidZonePressure
        );

    const label solidPatchStart = 
        stressMesh.boundaryMesh()[solidPatchID].start();

    forAll(solidPatchTraction, i)
    {
        solidPatchTraction[i] =
            solidZoneTraction
            [
                stressMesh.faceZones()[solidZoneID]
               .whichFace(solidPatchStart + i)
            ];
    }

    forAll(solidPatchPressure, i)
    {
        solidPatchPressure[i] =
            solidZonePressure
            [
                stressMesh.faceZones()[solidZoneID]
               .whichFace(solidPatchStart + i)
            ];
    }

    if (fsi)
    {
        tForce.traction() = solidPatchTraction;
        tForce.pressure() = solidPatchPressure;
    }

    vector totalTractionForce =
        gSum
        (
            solidPatchTraction
           *stressMesh.magSf().boundaryField()[solidPatchID]
        );

    Info << "Total traction force = "
        << totalTractionForce << endl;
}
