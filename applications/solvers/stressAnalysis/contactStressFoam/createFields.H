    mesh.weights();

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

    Info<< "Reading material properties\n" << endl;

    IOdictionary mechanicalProperties
    (
        IOobject
        (
            "mechanicalProperties",
            "constant",
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar rho
    (
        mechanicalProperties.lookup("rho")
    );

    dimensionedScalar rhoE
    (
        mechanicalProperties.lookup("E")
    );

    dimensionedScalar nu
    (
        mechanicalProperties.lookup("nu")
    );

    Info<< "Normalising E : E/rho\n" << endl;
    dimensionedScalar E = rhoE/rho;

    Info<< "Reading boundary conditions\n" << endl;

    IOdictionary boundaryConditions
    (
        IOobject
        (
            "boundaryConditions",
            "constant",
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    Info<< "Reading patchPressure\n" << endl;
    scalarField rhoPatchPressure
    (
        boundaryConditions.lookup("patchPressure")
    );

    Info<< "Normalising patchPressure : patchPressure/rho\n" << endl;
    scalarField patchPressure = rhoPatchPressure/rho.value();

    Info<< "Reading patchTraction\n" << endl;
    vectorField rhoPatchTraction
    (
        boundaryConditions.lookup("patchTraction")
    );

    Info<< "Normalising patchTraction : patchTraction/rho\n" << endl;
    vectorField patchTraction = rhoPatchTraction/rho.value();

    Info<< "Calculating Lame's coefficients\n" << endl;

    dimensionedScalar mu = E/(2.0*(1.0 + nu));
    dimensionedScalar lambda = nu*E/((1.0 + nu)*(1.0 - 2.0*nu));
    dimensionedScalar threeK = E/(1.0 - 2.0*nu);

    Switch planeStress
    (
        mechanicalProperties.lookup("planeStress")
    );

    if (planeStress)
    {
        Info<< "Plane Stress\n" << endl;

        //- change lambda and threeK for plane stress
        lambda = nu*E/((1.0 + nu)*(1.0 - nu));
    }
    else
    {
        Info<< "Plane Strain\n" << endl;
    }

    Info<< "mu = " << mu.value() << " Pa/rho\n";
    Info<< "lambda = " << lambda.value() << " Pa/rho\n";

#   include "calculatePatchTraction.H"

    Info<< "Reading contact patch IDs etc.\n" << endl;
    const dictionary& contactPatchDict
    (
        mesh.fvSolution::subDict("contactPatch")
    );
    label gradPatch(readLabel(contactPatchDict.lookup("gradPatchID")));
    label dirPatch(readLabel(contactPatchDict.lookup("dirPatchID")));
    scalar urf(readScalar(contactPatchDict.lookup("contactUnderrelaxation")));
    scalar touchTolerance
    (
        readScalar(contactPatchDict.lookup("touchTolerance"))
    );

    Info << "Creating contact pairs" << endl;
    contactPatchPair contactPair(mesh, gradPatch, dirPatch, touchTolerance);

    contactPatchPair reversePair(mesh, dirPatch, gradPatch, touchTolerance);

    volScalarField contactArea
    (
        IOobject
        (
            "contactArea",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(0)
    );

