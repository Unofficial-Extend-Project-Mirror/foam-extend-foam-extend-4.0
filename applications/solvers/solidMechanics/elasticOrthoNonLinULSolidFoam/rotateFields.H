//--------------------------------------------------//
//- rotate fields
//--------------------------------------------------//
{
    Info << "Rotating fields" << endl;

    volTensorField F = I + gradDU;

    U += DU;

    epsilon += DEpsilon;

    sigma += DSigma;

    volTensorField Finv = inv(F);

    volScalarField J = det(F);
    if(min(J.internalField()) < 0)
    {
        FatalErrorIn(args.executable())
            << "Negative Jacobian - a cell volume has become negative!"
            << exit(FatalError);
    }

    rho = rho/J;

    n = mesh.Sf()/mesh.magSf();

    //- rotate strain
    //epsilon = symm(Finv & epsilon & Finv.T());
    epsilon = transform(Finv, epsilon);

    //- rotate stress
    //sigma  = 1/J * symm(F.T() & sigma & F);
    sigma = (1/J) * transform(F.T(), sigma);

    //- rotate elastic constitutive tensor
    C = transform(F.T(), C);

    // - update implicit stiffness tensor
    forAll(K, celli)
    {
        K[celli].xx() = C[celli].xxxx();
        K[celli].yy() = C[celli].yyyy();
        K[celli].zz() = C[celli].zzzz();
    }
    K.correctBoundaryConditions();
}
