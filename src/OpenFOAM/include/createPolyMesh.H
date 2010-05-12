    Info<< "Create polyMesh for time = "
        << runTime.timeName() << nl << endl;

    polyMesh mesh
    (
        IOobject
        (
            polyMesh::defaultRegion,
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ
        )
    );
