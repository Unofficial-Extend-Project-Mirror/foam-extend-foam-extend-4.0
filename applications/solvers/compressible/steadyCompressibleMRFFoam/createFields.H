    Info<< "Reading thermophysical properties\n" << endl;

    autoPtr<basicPsiThermo> pThermo
    (
        basicPsiThermo::New(mesh)
    );
    basicPsiThermo& thermo = pThermo();

    volScalarField& p = thermo.p();
    volScalarField& h = thermo.h();
    const volScalarField& T = thermo.T();
    volScalarField psis("psi", thermo.psi()/thermo.Cp()*thermo.Cv());
    psis.oldTime();

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

    // Store kinetic energy for total enthalpy equation
    volScalarField K
    (
        "K",
        0.5*(magSqr(U))
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
            thermo
        )
    );

    // Create MRF zones
    MRFZones mrfZones(mesh);
    mrfZones.correctBoundaryVelocity(U);

    // Create relative velocity
    volVectorField Urel
    (
        IOobject
        (
            "Urel",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        U
    );
    mrfZones.relativeVelocity(Urel);

    // Create rotational velocity (= omega x r)
    volVectorField Urot
    (
        "Urot",
        U - Urel
    );

    // Create rothalpy, in two steps to preserve boundary conditions
    volScalarField i
    (
        IOobject
        (
            "i",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        h
    );
    i == h - 0.5*(magSqr(Urot) - magSqr(Urel));

    mesh.schemesDict().setFluxRequired(p.name());
