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

    Info<< "Reading solid thermal properties" << endl;
    thermalModel solidThermo(Tsolid);

    Info<< "Reading solid diffusivity k\n" << endl;
    volScalarField ksolid
    (
        IOobject
        (
            "k",
            runTime.timeName(),
            solidMesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        solidMesh
    );

    Info<< "Constructing rhoCp\n" << endl;
    volScalarField rhoCpsolid
    (
        IOobject
        (
            "rhoCpsolid",
            runTime.timeName(),
            solidMesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        solidThermo.rho()*solidThermo.C()
    );
    rhoCpsolid.oldTime();
