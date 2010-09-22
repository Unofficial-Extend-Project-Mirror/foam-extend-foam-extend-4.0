    scalar totalVolume = sum(mesh.V()).value();

    Info<< "Reading thermophysical properties\n" << endl;

    autoPtr<basicPsiThermo> pThermo
    (
        basicPsiThermo::New(mesh)
    );
    basicPsiThermo& thermo = pThermo();

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
        thermo.rho()
    );
    rho.oldTime();

    volScalarField& p = thermo.p();
    p.oldTime();

    const volScalarField& psi = thermo.psi();
    const volScalarField& T = thermo.T();
    volScalarField& h = thermo.h();


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
    autoPtr<compressible::turbulenceModel> turbulence
    (
        compressible::turbulenceModel::New
        (
            rho,
            U,
            phi,
            thermo
        )
    );

    Info<< "Creating field DpDt\n" << endl;
    volScalarField DpDt
    (
        IOobject
        (
            "DpDt",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::DDt(surfaceScalarField("phiU", phi/fvc::interpolate(rho)), p)
    );


    volScalarField rUA
    (
        IOobject
        (
            "rUA",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("rUA", dimensionSet(-1, 3, 1, 0, 0, 0, 0), 1),
        zeroGradientFvPatchScalarField::typeName
    );

    volVectorField H
    (
        IOobject
        (
            "H",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedVector("H", U.dimensions()/rUA.dimensions(), vector::zero),
        zeroGradientFvPatchVectorField::typeName
    );

    volScalarField divPhi
    (
        IOobject
        (
            "divPhi",
            runTime.timeName(),
            mesh
        ),
        fvc::div(phi)
    );


    surfaceScalarField meshFlux
    (
        IOobject
        (
            "meshFlux",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::interpolate(rho) * fvc::meshPhi(U)
    );
