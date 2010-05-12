    Info<< "Reading field Tsolid\n" << endl;
    volScalarField Tsolid
    (
        IOobject
        (
            "T",
            runTime.timeName(),
            solidMesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        solidMesh
    );

    Info<< "Reading solid diffusivity DT\n" << endl;

    volScalarField DTsolid
    (
        IOobject
        (
            "DT",
            runTime.timeName(),
            solidMesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        solidMesh
    );

    // On creation of both diffusivities, couple the sides
    DT.correctBoundaryConditions();
    DTsolid.correctBoundaryConditions();
