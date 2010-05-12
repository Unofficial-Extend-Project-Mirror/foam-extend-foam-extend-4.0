    fvMesh stressMesh
    (
        IOobject
        (
            "solid",
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ
        )
    );

    pointMesh pStressMesh(stressMesh);

    volPointInterpolation cpi
    (
        stressMesh,
        pStressMesh
    );

