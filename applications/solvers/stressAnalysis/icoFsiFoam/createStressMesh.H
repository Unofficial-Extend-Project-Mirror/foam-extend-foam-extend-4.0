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

    const volPointInterpolation& cpi = volPointInterpolation::New(stressMesh);

