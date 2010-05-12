    Info<< "Reading transportProperties\n" << endl;

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar nu
    (
        transportProperties.lookup("nu")
    );


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


    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, mesh.solutionDict().subDict("PISO"), pRefCell, pRefValue);

    scalar totalVolume = sum(mesh.V()).value();

    volScalarField rUA
    (
        IOobject
        (
            "rUA",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("rUA", dimTime, 1),
        zeroGradientFvPatchScalarField::typeName
    );


    {
        const dictionary& piso = mesh.solutionDict().subDict("PISO");

        bool momentumPredictor = true;
        if (piso.found("momentumPredictor"))
        {
            momentumPredictor = Switch(piso.lookup("momentumPredictor"));
        }

        bool momentumPredictorSave = momentumPredictor;
        momentumPredictor = false;
#       include "UEqn.H"
        momentumPredictor = momentumPredictorSave;
        
        rUA = 1.0/UEqn.A();
    }
