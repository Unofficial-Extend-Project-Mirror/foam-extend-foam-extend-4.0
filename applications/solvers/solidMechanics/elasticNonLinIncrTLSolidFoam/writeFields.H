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

    Info<< "Max epsilonEq = " << max(epsilonEq).value()
        << endl;

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

    //- Calculate Cauchy stress
    volTensorField F = I + gradU;
    volScalarField J = det(F);

    //- update density
    rho = rho/J;

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
        (1/J) * symm(F.T() & sigma & F)
    );

    //- Cauchy von Mises stress
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

    Info<< "Max sigmaCauchyEq = " << max(sigmaCauchyEq).value()
        << endl;

    volTensorField Finv = inv(F);
    volSymmTensorField epsilonAlmansi
    (
        IOobject
        (
            "epsilonAlmansi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
       ),
       symm(Finv & epsilon & Finv.T())
    );

    //- boundary traction
    volVectorField traction
    (
        IOobject
        (
            "tractionCauchy",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimForce/dimArea, vector::zero)
    );
    forAll(traction.boundaryField(), patchi)
    {
        tensorField Fpatch = I + gradU.boundaryField()[patchi];

        traction.boundaryField()[patchi] =
            n.boundaryField()[patchi] & (sigma.boundaryField()[patchi] & Fpatch);
    }

    runTime.write();
}
