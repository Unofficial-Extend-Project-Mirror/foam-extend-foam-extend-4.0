    Info<< "Reading field DU\n" << endl;
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

    volTensorField gradU
    (
        IOobject
        (
            "grad(U)",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedTensor("zero", dimless, tensor::zero)
    );

    //- increment of Green finite strain tensor
    volSymmTensorField DEpsilon
    (
        IOobject
        (
            "DEpsilonGreen",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    );

    //- Green strain tensor
    volSymmTensorField epsilon
    (
        IOobject
        (
            "epsilonGreen",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    );


    //- increment of 2nd Piola-Kirchhoff stress tensor
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

    //- 2nd Piola-Kirchhoff stress tensor
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

    constitutiveModel rheology(sigma, DU);

    volScalarField rho = rheology.rho();
    volScalarField mu = rheology.mu();
    volScalarField lambda = rheology.lambda();
    surfaceScalarField muf = fvc::interpolate(mu, "mu");
    surfaceScalarField lambdaf = fvc::interpolate(lambda, "lambda");

    surfaceVectorField n = mesh.Sf()/mesh.magSf();
