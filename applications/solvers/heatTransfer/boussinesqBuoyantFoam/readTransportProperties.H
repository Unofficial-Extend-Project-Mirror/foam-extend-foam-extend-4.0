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

    dimensionedScalar rho0
    (
        transportProperties.lookup("rho0")
    );

    dimensionedScalar T0
    (
        transportProperties.lookup("T0")
    );

    dimensionedScalar nu
    (
        dimensionedScalar(transportProperties.lookup("mu"))/rho0
    );

    dimensionedScalar DT
    (
        dimensionedScalar(transportProperties.lookup("k"))/
        (rho0*dimensionedScalar(transportProperties.lookup("Cp")))
    );

    dimensionedScalar beta
    (
        transportProperties.lookup("beta")
    );


