    Info<< "Create mesh for time = "
        << runTime.timeName() << nl << endl;

Info << "f1" << endl;

    fvMesh mesh
    (
        IOobject
        (
            fvMesh::defaultRegion,
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ
        )
    );

Info << "f2" << endl;

