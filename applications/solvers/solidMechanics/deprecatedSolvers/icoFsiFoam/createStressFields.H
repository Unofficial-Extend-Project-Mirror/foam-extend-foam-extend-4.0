    Info<< "Reading field U\n" << endl;
    volVectorField Usolid
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            stressMesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        stressMesh
    );
