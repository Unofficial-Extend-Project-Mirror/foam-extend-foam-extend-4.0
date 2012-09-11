vectorField fluidPatchPointsDispl
(
    mesh.boundaryMesh()[fluidPatchID].nPoints(),
    vector::zero
);

vectorField fluidPatchPointsDisplOld
(
    mesh.boundaryMesh()[fluidPatchID].nPoints(),
    vector::zero
);

vectorField solidPatchPointsDispl
(
    mesh.boundaryMesh()[fluidPatchID].nPoints(),
    vector::zero
);

vectorField fsiResidual
(
    mesh.boundaryMesh()[fluidPatchID].nPoints(),
    vector::zero
);

vectorField fsiResidualOld
(
    mesh.boundaryMesh()[fluidPatchID].nPoints(),
    vector::zero
);

vectorField solidPatchTraction
(
    stressMesh.boundary()[solidPatchID].size(), 
    vector::zero
);

scalarField solidPatchPressure
(
    stressMesh.boundary()[solidPatchID].size(), 
    0.0
);

scalar initialFsiResidualNorm = 0;
scalar fsiResidualNorm = 0;


