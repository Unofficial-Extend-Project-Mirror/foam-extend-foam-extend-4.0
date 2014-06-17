
//    Info<< "Reading field fvmU\n" << endl;

    pointMesh pMesh(mesh);

    pointVectorField fvmU
    (
        IOobject
        (
            "fvmU",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        pMesh,
        dimensionedVector("fvmU", dimLength, vector::zero),
        fvPatchVectorField::calculatedType()
    );
