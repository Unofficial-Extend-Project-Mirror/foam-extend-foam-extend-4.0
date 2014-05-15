    if (runTime.outputTime())
    {
        // Displacement gradient
        tetPointTensorField gradU = tetFec::grad(U);

        // Stress tensor
        tetPointSymmTensorField sigma =
            rho*(2.0*mu*symm(gradU) + lambda*I*tr(gradU));


        // Create pointMesh for field post-processing
        const pointMesh& pMesh = pointMesh::New(mesh);

        // U
        pointVectorField Up
        (
            IOobject
            (
                "Up",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            pMesh,
            U.dimensions()
        );

        Up.internalField() = vectorField::subField
        (
            U.internalField(),
            pMesh.size()
        );

        // sigmaEq
        pointScalarField sigmaEq
        (
            IOobject
            (
                "sigmaEq",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            pMesh,
            sigma.dimensions()
        );

        sigmaEq.internalField() = scalarField::subField
        (
            sqrt((3.0/2.0)*magSqr(dev(sigma.internalField())))(),
            pMesh.size()
        );

        // sigmaXX
        pointScalarField sigmaXX
        (
            IOobject
            (
                "sigmaXX",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            pMesh,
            sigma.dimensions()
        );

        sigmaXX.internalField() = scalarField::subField
        (
            sigma.component(symmTensor::XX)().internalField(),
            pMesh.size()
        );

        // sigmaYY
        pointScalarField sigmaYY
        (
            IOobject
            (
                "sigmaYY",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            pMesh,
            sigma.dimensions()
        );

        sigmaYY.internalField() = scalarField::subField
        (
            sigma.component(symmTensor::YY)().internalField(),
            pMesh.size()
        );

        // sigmaXY
        pointScalarField sigmaXY
        (
            IOobject
            (
                "sigmaXY",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            pMesh,
            sigma.dimensions()
        );

        sigmaXY.internalField() = scalarField::subField
        (
            sigma.component(symmTensor::XY)().internalField(),
            pMesh.size()
        );

        runTime.write();
    }
