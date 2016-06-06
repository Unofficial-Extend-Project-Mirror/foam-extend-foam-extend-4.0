    fvVectorMatrix UEqn
    (
        betav*fvm::ddt(rho, U)
      + fvm::div(phi, U)
      + turbulence->divDevRhoReff()
     ==
        betav*rho*g
    );

    volSymmTensorField invA = inv(I*UEqn.A() + drag->Dcu());

    if (piso.momentumPredictor())
    {
        U = invA & (UEqn.H() - betav*fvc::grad(p));
        U.correctBoundaryConditions();
    }
