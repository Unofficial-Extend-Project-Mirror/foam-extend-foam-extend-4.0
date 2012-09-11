{
    const vectorField& solidPatchDisplacement = 
        DU.boundaryField()[solidPatchID];

    vectorField solidZoneDisplacement
        (
            stressMesh.faceZones()[solidZoneID]().size(), 
            vector::zero
        );

    const label solidPatchStart = 
        stressMesh.boundaryMesh()[solidPatchID].start();

    forAll(solidPatchDisplacement, i)
    {
        solidZoneDisplacement
        [
            stressMesh.faceZones()[solidZoneID]
           .whichFace(solidPatchStart + i)
        ] = 
            solidPatchDisplacement[i];
    }

    // Parallel data exchange: collect displacement field on all processors
    reduce(solidZoneDisplacement, sumOp<vectorField>());

    vectorField fluidZoneDisplacement =
        interpolatorSolidFluid.faceInterpolate
        (
            solidZoneDisplacement
        );

    vectorField fluidPatchDisplacement
        (
            mesh.boundary()[fluidPatchID].size(),
            vector::zero
        );

    const label fluidPatchStart = 
        mesh.boundaryMesh()[fluidPatchID].start();

    forAll(fluidPatchDisplacement, i)
    {
        fluidPatchDisplacement[i] =
            fluidZoneDisplacement
            [
                mesh.faceZones()[fluidZoneID].whichFace(fluidPatchStart + i)
            ];
    }

    primitivePatchInterpolation fluidPatchInterpolator
        (
            mesh.boundaryMesh()[fluidPatchID]
        );

    solidPatchPointsDispl =
        fluidPatchInterpolator.faceToPointInterpolate
        (
            fluidPatchDisplacement
        );

    fsiResidualOld = fsiResidual;

    fsiResidual = solidPatchPointsDispl - fluidPatchPointsDispl;

//     maxFsiResidual = 
//        gMax
//         (
//             mag(fsiResidual)
//            /(mag(solidPatchPointsDispl) + SMALL)
//         );

//     Info << "Maximal fsi residual: " << maxFsiResidual << endl;

    fsiResidualNorm = ::sqrt(gSum(magSqr(fsiResidual)));

    if (outerCorr == 1)
    {
        initialFsiResidualNorm = fsiResidualNorm;
    }

    fsiResidualNorm /= initialFsiResidualNorm + SMALL;

    Info << "Current fsi residual norm: " << fsiResidualNorm << endl;
}
