    if (runTime.outputTime())
    {
        // Displacement gradient
        tetPointTensorField gradU = tetFec::grad(U);

        // Stress tensor
        tetPointSymmTensorField sigma =
            rho*(2.0*mu*symm(gradU) + lambda*I*tr(gradU));

        // sigmaXX
        tetPointScalarField sigmaXX
        (
            IOobject
            (
                "sigmaXX",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            sigma.component(symmTensor::XX)
        );

        // sigmaYY
        tetPointScalarField sigmaYY
        (
            IOobject
            (
                "sigmaYY",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            sigma.component(symmTensor::YY)
        );

        // sigmaXY
        tetPointScalarField sigmaXY
        (
            IOobject
            (
                "sigmaXY",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            sigma.component(symmTensor::XY)
        );

        runTime.write();
    }
