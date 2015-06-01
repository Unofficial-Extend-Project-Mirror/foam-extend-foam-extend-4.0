    Info<< "Reading incremental displacement field DU\n" << endl;
    volVectorField DU
    (
        IOobject
        (
            "DU",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volTensorField gradDU = fvc::grad(DU);

    Info<< "Reading accumulated displacement field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimLength, vector::zero)
    );

    Info << "Reading accumulated strain field epsilon\n" << endl;
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

    // total plastic strain
    volSymmTensorField epsilonP
    (
        IOobject
        (
            "epsilonP",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    );

    volSymmTensorField DEpsilon
    (
        IOobject
        (
            "DEpsilon",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    );

    Info << "Reading accumulated stress field sigma\n" << endl;
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


    Info << "Reading incremental stress field DSigma\n" << endl;
    volSymmTensorField DSigma
    (
        IOobject
        (
            "DSigma",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimForce/dimArea, symmTensor::zero)
    );

    //- explicit terms in the momentum equation
    volVectorField divDSigmaExp
    (
        IOobject
        (
            "divDSigmaExp",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimensionSet(1, -2, -2, 0, 0, 0, 0), vector::zero)
    );

    volVectorField divDSigmaNonLinExp
    (
        IOobject
        (
            "divDSigmaNonLinExp",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimensionSet(1,-2,-2,0,0,0,0), vector::zero)
    );

    constitutiveModel rheology(sigma, DU);

    volScalarField rho = rheology.rho();
    volScalarField mu = rheology.mu();
    volScalarField lambda = rheology.lambda();
    surfaceScalarField muf = fvc::interpolate(rheology.mu());
    surfaceScalarField lambdaf = fvc::interpolate(rheology.lambda());

    surfaceVectorField n = mesh.Sf()/mesh.magSf();

    // plastic strain increment
    const volSymmTensorField& DEpsilonP = rheology.DEpsilonP();

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
    scalar aitkenTheta = 0.1;

//    volVectorField resid
//    (
//        IOobject
//        (
//            "resid",
//            runTime.timeName(),
//            mesh,
//            IOobject::NO_READ,
//            IOobject::AUTO_WRITE
//        ),
//        mesh,
//        dimensionedVector("zero", dimless, vector::zero)
//    );
