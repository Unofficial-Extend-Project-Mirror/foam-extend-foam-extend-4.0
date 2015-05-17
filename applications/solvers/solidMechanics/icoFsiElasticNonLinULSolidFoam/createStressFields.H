    Info<< "Reading incremental displacement field DU\n" << endl;

    volVectorField DU
    (
        IOobject
        (
            "DU",
            runTime.timeName(),
            stressMesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        stressMesh
    );

    volTensorField gradDU = fvc::grad(DU);

    volVectorField Usolid
    (
        IOobject
        (
            "Usolid",
            runTime.timeName(),
            stressMesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        DU
    );


    Info<< "Reading incremental displacement field DV\n" << endl;
    volVectorField DV
    (
        IOobject
        (
            "DV",
            runTime.timeName(),
            stressMesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::ddt(DU)
    );


    Info<< "Reading accumulated velocity field V\n" << endl;
    volVectorField Vs
    (
        IOobject
        (
            "Vs",
            runTime.timeName(),
            stressMesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        stressMesh,
        dimensionedVector("zero", dimVelocity, vector::zero)
    );


    Info << "Reading accumulated stress field sigma\n" << endl;
    volSymmTensorField sigma
    (
        IOobject
        (
            "sigma",
            runTime.timeName(),
            stressMesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        stressMesh,
        dimensionedSymmTensor("zero", dimForce/dimArea, symmTensor::zero)
    );


    Info << "Reading incremental stress field DSigma\n" << endl;
    volSymmTensorField DSigma
    (
        IOobject
        (
            "DSigma",
            runTime.timeName(),
            stressMesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        stressMesh,
        dimensionedSymmTensor("zero", dimForce/dimArea, symmTensor::zero)
    );

    constitutiveModel rheology(sigma, DU);
    volScalarField rho = rheology.rho();
    volScalarField mu = rheology.mu();
    volScalarField lambda = rheology.lambda();

    volTensorField F = I + gradDU.T();
    volTensorField DF = F - I;
    volScalarField J = det(F);


    word solidDdtScheme
    (
        stressMesh.schemesDict().ddtScheme("ddt(" + DU.name() +')')
    );

//     if
//     (
//         solidDdtScheme != fv::EulerFixedMeshDdtScheme<vector>::typeName
//      && solidDdtScheme != fv::backwardFixedMeshDdtScheme<vector>::typeName
//     )
//     {
//         FatalErrorIn(args.executable())
//             << "Selected temporal differencing scheme: " << solidDdtScheme
//             << ", instead: "
//             << fv::EulerFixedMeshDdtScheme<vector>::typeName
//             << " of "
//             << fv::backwardFixedMeshDdtScheme<vector>::typeName
//             << exit(FatalError);
//     }


//     IOdictionary rheologyProp
//     (
//         IOobject
//         (
//             "rheologyProperties",
//             runTime.constant(),
//             stressMesh,
//             IOobject::MUST_READ,
//             IOobject::NO_WRITE
//         )
//     );


// dimensionedVector fb(rheologyProp.lookup("bodyForce"));
