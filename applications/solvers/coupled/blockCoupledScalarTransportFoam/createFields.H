    Info<< "Reading field T\n" << endl;
    volScalarField T
    (
        IOobject
        (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    Info<< "Reading field Ts\n" << endl;
    volScalarField Ts
    (
        IOobject
        (
            "Ts",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


// Removed comparison with the reference field.  HJ, 17/Jun/2010
#if (0)
    Info<< "Reading field Tref\n" << endl;
    volScalarField Tref
    (
        IOobject
        (
            "Tref",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    Info<< "Reading field Tsref\n" << endl;
    volScalarField Tsref
    (
        IOobject
        (
            "Tsref",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
#endif

    Info<< "Reading field U\n" << endl;

    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


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


    Info<< "Reading diffusivity D\n" << endl;

    dimensionedScalar DT
    (
        transportProperties.lookup("DT")
    );

    dimensionedScalar DTs
    (
        transportProperties.lookup("DTs")
    );

    dimensionedScalar alpha
    (
        transportProperties.lookup("alpha")
    );

#   include "createPhi.H"
