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

#   include "createPhi.H"

    // Set reference cell taking into account whether implicit or explicit
    // overset aglorithm is used for field p.
    const oversetMesh& om = oversetMesh::New(mesh);
    label pRefCell = 0;
    scalar pRefValue = 0.0;
    om.setRefCell(p, pimple.dict(), pRefCell, pRefValue);

    mesh.schemesDict().setFluxRequired(p.name());
    if (correctPhi)
    {
        mesh.schemesDict().setFluxRequired("pcorr");
    }

    singlePhaseTransportModel laminarTransport(U, phi);

    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(U, phi, laminarTransport)
    );

    Info<< "Reading field rAU if present\n" << endl;
    volScalarField rAU
    (
        IOobject
        (
            "rAU",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        runTime.deltaT(),
        zeroGradientFvPatchScalarField::typeName
    );

#   include "writeOversetMasks.H"
