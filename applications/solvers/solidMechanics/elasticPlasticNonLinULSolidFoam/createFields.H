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

    Info<< "Creating field U\n" << endl;
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

     volVectorField divDSigmaExp
	(
	 IOobject
	 (
	  "divDSigmaExp",
	  runTime.timeName(),
	  mesh,
	  IOobject::NO_READ,
	  IOobject::AUTO_WRITE
	  ),
	 mesh,
 	 dimensionedVector("zero", dimensionSet(1, -2, -2, 0, 0, 0, 0), vector::zero)
	 );

     volVectorField divDSigmaLargeStrainExp
	(
	 IOobject
	 (
	  "divDSigmaLargeStrainExp",
	  runTime.timeName(),
	  mesh,
	  IOobject::NO_READ,
	  IOobject::AUTO_WRITE
	  ),
	 mesh,
 	 dimensionedVector("zero", dimensionSet(1, -2, -2, 0, 0, 0, 0), vector::zero)
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

    volSymmTensorField epsilonE
    (
        IOobject
        (
            "epsilonE",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("zero", dimless, symmTensor::zero)
    );

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

    plasticityModel rheology(gradDU, epsilon, sigma);

    volScalarField mu = rheology.newMu();
    volScalarField lambda = rheology.newLambda();
    surfaceScalarField muf = fvc::interpolate(rheology.newMu());
    surfaceScalarField lambdaf = fvc::interpolate(rheology.newLambda());

    surfaceVectorField n = mesh.Sf()/mesh.magSf();

    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        rheology.rho()
     );

    volTensorField DF = gradDU.T();
    volTensorField F = I + DF;
    volScalarField J = det(F);
