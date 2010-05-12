    Info<< "Reading thermophysical properties\n" << endl;

    autoPtr<basicThermo> thermo
    (
        basicThermo::New(mesh)
    );

    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        thermo->rho()
    );

    volScalarField& p = thermo->p();
    volScalarField& h = thermo->h();
    const volScalarField& psi = thermo->psi();


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

    #include "compressibleCreatePhi.H"


    Info<< "Creating turbulence model\n" << endl;
    autoPtr<compressible::LESModel> turbulence
    (
        compressible::LESModel::New(rho, U, phi, thermo())
    );


    Info<< "Creating field dpdt\n" << endl;
    volScalarField dpdt = fvc::ddt(p);

    Info<< "Calculating field g.h\n" << endl;
    volScalarField gh("gh", g & mesh.C());

    dimensionedScalar pRef("pRef", p.dimensions(), thermo->lookup("pRef"));

    Info<< "Creating field pd\n" << endl;
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

    p = pd + rho*gh + pRef;
    thermo->correct();

    dimensionedScalar initialMass = fvc::domainIntegrate(rho);
