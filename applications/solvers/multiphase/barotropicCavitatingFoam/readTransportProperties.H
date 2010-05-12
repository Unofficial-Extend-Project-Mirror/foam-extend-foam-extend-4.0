    Info<< "Reading transportProperties\n" << endl;

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar nug
    (
        transportProperties.lookup("nug")
    );

    dimensionedScalar nul
    (
        transportProperties.lookup("nul")
    );
