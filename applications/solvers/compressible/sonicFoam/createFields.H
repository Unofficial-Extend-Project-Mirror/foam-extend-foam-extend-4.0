    Info<< "Reading thermophysical properties\n" << endl;

    autoPtr<basicPsiThermo> pThermo
    (
        basicPsiThermo::New(mesh)
    );
    basicPsiThermo& thermo = pThermo();

    volScalarField& p = thermo.p();
    volScalarField& e = thermo.e();

    volScalarField psis("psi", thermo.psi()/thermo.Cp()*thermo.Cv());
    psis.oldTime();

    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh
        ),
        thermo.rho()
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

    mesh.schemesDict().setFluxRequired(p.name());
