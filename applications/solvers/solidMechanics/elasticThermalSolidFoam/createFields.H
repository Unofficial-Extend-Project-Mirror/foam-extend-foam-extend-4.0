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

    volTensorField gradU
    (
        IOobject
        (
            "grad(U)",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedTensor("zero", dimless, tensor::zero)
    );

    surfaceVectorField snGradU
    (
        IOobject
        (
            "snGrad(U)",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimless, vector::zero)
    );

    volSymmTensorField epsilon
    (
        IOobject
        (
            "epsilon",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    );

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

    volSymmTensorField sigma
    (
        IOobject
        (
            "sigma",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimForce/dimArea, symmTensor::zero)
    );

    volVectorField divSigmaExp
    (
        IOobject
        (
            "divSigmaExp",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimensionSet(1,-2,-2,0,0,0,0), vector::zero)
    );

    surfaceVectorField n = mesh.Sf()/mesh.magSf();

    // mechanical properties
    constitutiveModel rheology(sigma, U);

    volScalarField rho = rheology.rho();

    volScalarField mu = rheology.mu();
    volScalarField lambda = rheology.lambda();
    surfaceScalarField muf = fvc::interpolate(mu, "mu");
    surfaceScalarField lambdaf = fvc::interpolate(lambda, "lambda");


    // thermal properties
    Info<< "Reading thermal model\n" << endl;
    thermalModel thermal(T);

    volScalarField C = thermal.C();
    volScalarField k
    (
        "DT",
        thermal.k()
    );
    volScalarField threeKalpha = rheology.threeK()*rho*thermal.alpha();
    surfaceScalarField threeKalphaf =
        fvc::interpolate(threeKalpha, "threeKalpha");
    volScalarField T0 = thermal.T0();
    volScalarField rhoC = rho*C;


    // for aitken relaxation
    volVectorField aitkenDelta
    (
        IOobject
        (
         "aitkenDelta",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE
         ),
        mesh,
        dimensionedVector("zero", dimLength, vector::zero)
    );

    // aitken relaxation factor
    scalar aitkenInitialRes = 1.0;
    scalar aitkenTheta = 0.01;
