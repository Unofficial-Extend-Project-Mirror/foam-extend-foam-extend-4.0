// Read velocity field for turbulence
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

// Read temperature field
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// Create the flux for turbulence
#include "createPhi.H"

// Read transport properties and all necessary coefficients
#include "readTransportProperties.H"

// Create generic turbulence model (RAS/LES)
autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New
    (
        U,
        phi,
        laminarTransport
    )
);
