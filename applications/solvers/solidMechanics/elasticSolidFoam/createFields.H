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
        dimensionedVector("zero", dimForce/dimVolume, vector::zero)
    );

    constitutiveModel rheology(sigma, U);

    volScalarField rho = rheology.rho();
    volScalarField mu = rheology.mu();
    volScalarField lambda = rheology.lambda();
    surfaceScalarField muf = rheology.muf();
    surfaceScalarField lambdaf = rheology.lambdaf();

    surfaceVectorField n = mesh.Sf()/mesh.magSf();

// check if there are contact boundaries to know
// if contact fields should be written for visualisation
// note: this is just for visualisation
// bool contactBoundaries = false;
// forAll(U.boundaryField(), patchi)
// {
//   if(U.boundaryField()[patchi].type()
//      ==
//      solidContactFvPatchVectorField::typeName)
//     {
//       contactBoundaries = true;
//       break;
//     }
// }


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
