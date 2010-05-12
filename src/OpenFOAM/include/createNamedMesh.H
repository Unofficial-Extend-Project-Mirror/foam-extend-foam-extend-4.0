    word regionName;

    if (args.options().found("region"))
    {
        regionName = args.options()["region"];

        Info<< "Create mesh " << regionName << " for time = "
            << runTime.timeName() << nl << endl;
    }
    else
    {
        regionName = fvMesh::defaultRegion;

        Info<< "Create mesh for time = "
            << runTime.timeName() << nl << endl;
    }

    fvMesh mesh
    (
        IOobject
        (
            regionName,
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ
        )
    );
