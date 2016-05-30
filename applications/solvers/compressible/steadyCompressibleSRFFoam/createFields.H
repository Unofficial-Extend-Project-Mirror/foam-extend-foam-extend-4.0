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

    Info<< "\nReading field Urel\n" << endl;
    volVectorField Urel
    (
        IOobject
        (
            "Urel",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading/calculating face flux field phi\n" << endl;

    surfaceScalarField phi
    (
        IOobject
        (
            "phi",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        linearInterpolate(rho*Urel) & mesh.Sf()
    );

    Info<< "Creating turbulence model\n" << endl;
    autoPtr<compressible::RASModel> turbulence
    (
        compressible::RASModel::New
        (
            rho,
            Urel,
            phi,
            thermo
        )
    );

    Info<< "Creating SRF model\n" << endl;
    autoPtr<SRF::SRFModel> SRF
    (
        SRF::SRFModel::New(Urel)
    );

    // Create rotational velocity
    volVectorField Urot("Urot", SRF->U());

    // Create absolute velocity field
    volVectorField Uabs
    (
        IOobject
        (
            "Uabs",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        Urel + Urot
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
