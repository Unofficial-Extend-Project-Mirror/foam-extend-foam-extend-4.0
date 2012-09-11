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
    surfaceVectorField snGradU = fvc::snGrad(U);

    volVectorField V
    (
        IOobject
        (
            "V",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::ddt(U)
    );

   volTensorField gradV = fvc::ddt(gradU);
   surfaceVectorField snGradV = (snGradU - snGradU.oldTime())/runTime.deltaT();

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


   // read rheology properties
    rheologyModel rheology(sigma);

    volScalarField rho = rheology.rho();

    volScalarField mu = rheology.mu();
    volScalarField lambda = rheology.lambda();
    surfaceScalarField muf = fvc::interpolate(mu, "mu");
    surfaceScalarField lambdaf = fvc::interpolate(lambda, "lambda");

    surfaceVectorField n = mesh.Sf()/mesh.magSf();


   //- create contact problem
   contactProblem contact(U);
