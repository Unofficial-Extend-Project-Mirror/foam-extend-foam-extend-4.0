    Info<< "Reading field pd\n" << endl;
    volScalarField pd
    (
        IOobject
        (
            "pd",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

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

#   include "createPhi.H"

    multiphaseMixture mixture(U, phi);

    // Need to store rho for ddt(rho, U)
    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT
        ),
        mixture.rho()
    );
    rho.oldTime();


    Info<< "Calculating field g.h\n" << endl;
    surfaceScalarField ghf("gh", g & mesh.Cf());


    label pdRefCell = 0;
    scalar pdRefValue = 0.0;
    setRefCell(pd, mesh.solutionDict().subDict("PISO"), pdRefCell, pdRefValue);
