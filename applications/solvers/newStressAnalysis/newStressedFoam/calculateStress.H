    if (runTime.outputTime())
    {
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

        Info<< "Max sigmaEq = " << max(sigmaEq).value()
            << endl;

        volScalarField sigmaxx
        (
            IOobject
            (
                "sigmaxx",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            sigma.component(symmTensor::XX)
        );

        volScalarField sigmayy
        (
            IOobject
            (
                "sigmayy",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            sigma.component(symmTensor::YY)
        );

        volScalarField sigmazz
        (
            IOobject
            (
                "sigmazz",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            sigma.component(symmTensor::ZZ)
        );

        Info<< "Max sigmazz = " << max(sigmazz).value()
            << endl;

        volScalarField sigmaxy
        (
            IOobject
            (
                "sigmaxy",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            sigma.component(symmTensor::XY)
        );

        volScalarField sigmaxz
        (
            IOobject
            (
                "sigmaxz",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            sigma.component(symmTensor::XZ)
        );

        volScalarField sigmayz
        (
            IOobject
            (
                "sigmayz",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            sigma.component(symmTensor::YZ)
        );

        runTime.write();
    }
