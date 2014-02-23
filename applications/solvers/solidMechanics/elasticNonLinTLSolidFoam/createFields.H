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

    volTensorField gradU = fvc::grad(U);

    //- Green finite strain tensor
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

    //- Second Piola-Kirchhoff stress tensor
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

    constitutiveModel rheology(sigma, U);

    volScalarField rho = rheology.rho();
    volScalarField mu = rheology.mu();
    volScalarField lambda = rheology.lambda();
    surfaceScalarField muf = fvc::interpolate(mu, "mu");
    surfaceScalarField lambdaf = fvc::interpolate(lambda, "lambda");

    surfaceVectorField n = mesh.Sf()/mesh.magSf();
