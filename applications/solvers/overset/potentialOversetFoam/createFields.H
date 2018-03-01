    Info<< "Reading field p\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    // Do not reset p

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

    if (args.optionFound("resetU"))
    {
        U.internalField() = vector::zero;
    }

    surfaceScalarField phi
    (
        IOobject
        (
            "phi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fvc::interpolate(U) & mesh.Sf()
    );


    // Set reference cell taking into account whether implicit or explicit
    // overset algorithm is used for field p.
    label pRefCell = 0;
    scalar pRefValue = 0.0;
    om.setRefCell(p, simple.dict(), pRefCell, pRefValue);

    mesh.schemesDict().setFluxRequired(p.name());
