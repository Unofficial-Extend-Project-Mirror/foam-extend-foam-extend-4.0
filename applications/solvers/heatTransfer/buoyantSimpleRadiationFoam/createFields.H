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
    const volScalarField& T = thermo->T();


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

#   include "compressibleCreatePhi.H"


    Info<< "Creating turbulence model\n" << endl;
    autoPtr<compressible::RASModel> turbulence
    (
        compressible::RASModel::New
        (
            rho,
            U,
            phi,
            thermo()
        )
    );

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


    label pdRefCell = 0;
    scalar pdRefValue = 0.0;
    setRefCell
    (
        pd,
        mesh.solutionDict().subDict("SIMPLE"),
        pdRefCell,
        pdRefValue
    );


    Info<< "Creating radiation model\n" << endl;
    autoPtr<radiation::radiationModel> radiation
    (
        radiation::radiationModel::New(T)
    );


    dimensionedScalar initialMass = fvc::domainIntegrate(rho);
