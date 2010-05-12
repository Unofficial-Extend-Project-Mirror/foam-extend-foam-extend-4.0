    Info<< "Reading thermophysical properties\n" << endl;

    autoPtr<basicThermo> thermo
    (
        basicThermo::New(mesh)
    );

    volScalarField& p = thermo->p();
    volScalarField& h = thermo->h();
    const volScalarField& psi = thermo->psi();

    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        thermo->rho()
    );

    Info<< "\nReading field U\n" << endl;
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


    Info<< "Creating field DpDt\n" << endl;
    volScalarField DpDt =
        fvc::DDt(surfaceScalarField("phiU", phi/fvc::interpolate(rho)), p);

    pointMesh pMesh(mesh);
    volPointInterpolation vpi(mesh, pMesh);

    Info<< "Constructing thermoCloud1" << endl;
    basicThermoCloud thermoCloud1
    (
        "thermoCloud1",
        vpi,
        rho,
        U,
        g,
        thermo()
    );

    Info<< "Constructing kinematicCloud1" << endl;
    basicKinematicCloud kinematicCloud1
    (
        "kinematicCloud1",
        vpi,
        rho,
        U,
        thermo().mu(),
        g
    );

