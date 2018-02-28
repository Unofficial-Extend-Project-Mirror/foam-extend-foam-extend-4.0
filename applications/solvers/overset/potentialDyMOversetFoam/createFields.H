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

    surfaceScalarField phi
    (
        IOobject
        (
            "phi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fvc::interpolate(U) & mesh.Sf()
    );

    // Write overset div phi for visualisation purposes
    volScalarField oversetDivPhi
    (
        IOobject
        (
            "oversetDivPhi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fvc::div(phi)
    );

    // Set reference cell taking into account whether implicit or explicit
    // overset algorithm is used for field p.
    const oversetMesh& om = oversetMesh::New(mesh);
    label pRefCell = 0;
    scalar pRefValue = 0.0;
    om.setRefCell(p, piso.dict(), pRefCell, pRefValue);

    mesh.schemesDict().setFluxRequired(p.name());

    // Write overset masks
#   include "writeOversetMasks.H"
