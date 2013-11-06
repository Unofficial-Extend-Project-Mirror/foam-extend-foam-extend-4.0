{
    // Setting pressure on solid patch
    Info << "Setting pressure" << endl;

    scalarField solidPatchPressure =
        interpolatorFluidSolid.faceInterpolate
        (
            p.boundaryField()[fluidPatchID]
        );

    solidPatchPressure *= rhoFluid.value();
   
    tForce.pressure() = solidPatchPressure;


    vector totalPressureForce =
        sum
        (
            p.boundaryField()[fluidPatchID]*
            mesh.Sf().boundaryField()[fluidPatchID]
        );

     
    Info << "Total pressure force = " << totalPressureForce << endl;
}
