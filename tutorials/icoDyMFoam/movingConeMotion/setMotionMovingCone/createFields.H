    tetPolyMesh tetMesh(mesh);

    Info<< "Creating field motionGamma\n" << endl;
    elementScalarField motionGamma
    (
        IOobject
        (
            "motionGamma",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        tetMesh,
        dimensionedScalar("1.0", dimless, 1.0)
    );
