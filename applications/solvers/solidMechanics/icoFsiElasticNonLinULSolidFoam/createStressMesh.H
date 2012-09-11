// Info << " Reading stress mesh" << endl;

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

