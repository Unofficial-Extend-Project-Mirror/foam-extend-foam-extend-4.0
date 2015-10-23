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

    Info<< "Max epsilonPEq = " << max(epsilonPEq).value() << endl;

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

    volScalarField sigmaHyd
    (
        IOobject
        (
            "sigmaHyd",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        (
            sigma.component(symmTensor::XX)
          + sigma.component(symmTensor::YY)
          + sigma.component(symmTensor::ZZ)
        )/3
    );
    Info<< "Max sigmaHyd = " << max(sigmaHyd).value() << endl;

    runTime.write();
}
