    // Create volume-to surface mapping object
    volSurfaceMapping vsm(aMesh);

    volScalarField Cvf
    (
        IOobject
        (
            "Cvf",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimless/dimLength, 0)
    );

    vsm.mapToVolume(Cs, Cvf.boundaryField());
    Cvf.write();

    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimVelocity, vector::zero)
    );

    vsm.mapToVolume(Us, U.boundaryField());
    U.write();
