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

    Info<< "Creating phaseChangeTwoPhaseMixture\n" << endl;
    autoPtr<phaseChangeTwoPhaseMixture> twoPhaseProperties =
        phaseChangeTwoPhaseMixture::New(U, phi, "gamma");

    const dimensionedScalar& rho1 = twoPhaseProperties->rho1();
    const dimensionedScalar& rho2 = twoPhaseProperties->rho2();
    const dimensionedScalar& pSat = twoPhaseProperties->pSat();

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
        gamma*rho1 + (scalar(1) - gamma)*rho2,
        gamma.boundaryField().types()
    );
    rho.oldTime();


    label pdRefCell = 0;
    scalar pdRefValue = 0.0;
    setRefCell(pd, mesh.solutionDict().subDict("PISO"), pdRefCell, pdRefValue);

    Info<< "Calculating field g.h" << endl;
    volScalarField gh("gh", g & mesh.C());
    surfaceScalarField ghf("ghf", g & mesh.Cf());

    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        pd + rho*gh
    );


    // Construct interface from gamma distribution
    interfaceProperties interface(gamma, U, twoPhaseProperties());

    // Construct LES model
    autoPtr<incompressible::LESModel> turbulence
    (
        incompressible::LESModel::New(U, phi, twoPhaseProperties())
    );
