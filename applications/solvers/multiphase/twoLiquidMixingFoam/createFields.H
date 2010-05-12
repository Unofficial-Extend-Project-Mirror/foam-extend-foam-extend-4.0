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

    Info<< "Reading field gamma\n" << endl;
    volScalarField gamma
    (
        IOobject
        (
            "gamma",
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

    Info<< "Reading transportProperties\n" << endl;
    twoPhaseMixture twoPhaseProperties(U, phi, "gamma");
    
    const dimensionedScalar& rho1 = twoPhaseProperties.rho1();
    const dimensionedScalar& rho2 = twoPhaseProperties.rho2();

    dimensionedScalar Dab(twoPhaseProperties.lookup("Dab"));

    // Need to store rho for ddt(rho, U)
    volScalarField rho("rho", gamma*rho1 + (scalar(1) - gamma)*rho2);
    rho.oldTime();


    // Mass flux
    // Initialisation does not matter because rhoPhi is reset after the
    // gamma solution before it is used in the U equation.
    surfaceScalarField rhoPhi
    (
        IOobject
        (
            "rho*phi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rho1*phi
    );


    Info<< "Calculating field g.h\n" << endl;
    surfaceScalarField ghf("gh", g & mesh.Cf());


    label pdRefCell = 0;
    scalar pdRefValue = 0.0;
    setRefCell(pd, mesh.solutionDict().subDict("PISO"), pdRefCell, pdRefValue);
