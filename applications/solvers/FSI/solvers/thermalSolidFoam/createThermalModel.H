    Info<< "Reading field T\n" << endl;
    volScalarField T
    (
        IOobject
        (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    T.oldTime();


    // Thermal model
    Info<< "Reading thermal model\n" << endl;
    thermalModel thermal(T);

    volScalarField C = thermal.C();
    volScalarField k = thermal.k();
    volScalarField T0 = thermal.T0();

    volScalarField rho = stress().rheology().rho();
    volScalarField rhoC = rho*C;
