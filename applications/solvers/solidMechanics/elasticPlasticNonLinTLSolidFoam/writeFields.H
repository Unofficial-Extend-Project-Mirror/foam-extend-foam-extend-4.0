if (runTime.outputTime())
{
    volScalarField epsilonEq
    (
        IOobject
        (
            "epsilonEq",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        sqrt((2.0/3.0)*magSqr(dev(epsilon)))
    );
    Info<< "Max epsilonEq = " << max(epsilonEq).value() << endl;

    volScalarField epsilonPEq
    (
        IOobject
        (
            "epsilonPEq",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        sqrt((2.0/3.0)*magSqr(dev(epsilonP)))
    );
    Info<< "Max epsilonPEq = " << max(epsilonPEq).value()<< endl;

    volScalarField sigmaEq
    (
        IOobject
        (
            "sigmaEq",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        sqrt((3.0/2.0)*magSqr(dev(sigma)))
    );

    Info<< "Max sigmaEq = " << max(sigmaEq).value() << endl;

    // Deformation gradient
    volTensorField F = I + gradU;
    volScalarField J = det(F);

    // Calculate Cauchy stress
    volSymmTensorField sigmaCauchy
    (
        IOobject
        (
            "sigmaCauchy",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        (1/J)*symm(F.T() & sigma & F)
    );

    // Cauchy von Mises stress
    volScalarField sigmaCauchyEq
    (
        IOobject
        (
            "sigmaCauchyEq",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        sqrt((3.0/2.0)*magSqr(dev(sigmaCauchy)))
    );

    Info<< "Max sigmaCauchyEq = " << max(sigmaCauchyEq).value() << endl;

//     volTensorField Finv = inv(F);
//     volSymmTensorField epsilonAlmansi
//     (
//         IOobject
//         (
//             "epsilonAlmansi",
//             runTime.timeName(),
//             mesh,
//             IOobject::NO_READ,
//             IOobject::AUTO_WRITE
//         ),
//         symm(Finv & epsilon & Finv.T())
//     );

    runTime.write();
}
