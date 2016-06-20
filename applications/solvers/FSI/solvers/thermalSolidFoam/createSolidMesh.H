Info << " Reading stress mesh" << endl;

fvMesh mesh
(
    IOobject
    (
        Foam::fvMesh::defaultRegion,
        runTime.timeName(),
        runTime,
        IOobject::MUST_READ
    )
);
