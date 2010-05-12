    word regionName;

    if (args.options().found("region"))
    {
        regionName = args.options()["region"];

        Info<< "Create polyMesh " << regionName << " for time = "
            << runTime.timeName() << nl << endl;
    }
    else
    {
        regionName = polyMesh::defaultRegion;

        Info<< "Create polyMesh for time = "
            << runTime.timeName() << nl << endl;
    }

    polyMesh mesh
    (
        IOobject
        (
            regionName,
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ
        )
    );
