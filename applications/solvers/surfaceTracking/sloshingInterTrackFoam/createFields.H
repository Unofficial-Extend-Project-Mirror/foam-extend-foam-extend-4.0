Info<< "Reading field p\n" << endl;

volScalarField p
(
    IOobject
    (
        "p",
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


#include "createPhi.H"


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
    mesh,
    dimensionedScalar("0", dimMass/dimVolume, 0)
);

freeSurface interface(mesh, rho, U, p, phi);

rho = interface.fluidIndicator()*
    (
        interface.rhoFluidA() 
      - interface.rhoFluidB()
    ) 
  + interface.rhoFluidB();


volScalarField mu
(
    IOobject
    (
        "mu",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    interface.fluidIndicator()*
    (
        interface.muFluidA() 
      - interface.muFluidB()
    ) 
  + interface.muFluidB()
    );


surfaceScalarField phiNet = fvc::interpolate(rho)*(phi - fvc::meshPhi(rho,U));
phiNet.rename("phiNet");
