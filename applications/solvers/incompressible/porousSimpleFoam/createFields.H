    Info << "Reading field p\n" << endl;
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

    Info << "Reading field U\n" << endl;
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


    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, simple.dict(), pRefCell, pRefValue);
    mesh.schemesDict().setFluxRequired(p.name());

    singlePhaseTransportModel laminarTransport(U, phi);

    autoPtr<incompressible::RASModel> turbulence
    (
        incompressible::RASModel::New(U, phi, laminarTransport)
    );


    porousZones pZones(mesh);
    Switch pressureImplicitPorosity(false);

    int nUCorr = 0;
    if (pZones.size())
    {
        // nUCorrectors for pressureImplicitPorosity
        if (simple.dict().found("nUCorrectors"))
        {
            nUCorr = readInt(simple.dict().lookup("nUCorrectors"));
        }

        if (nUCorr > 0)
        {
            pressureImplicitPorosity = true;
        }
    }
