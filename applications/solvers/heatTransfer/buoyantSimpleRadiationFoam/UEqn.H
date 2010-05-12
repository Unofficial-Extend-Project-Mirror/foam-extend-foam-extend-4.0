    // Solve the Momentum equation

    tmp<fvVectorMatrix> UEqn
    (
        fvm::div(phi, U)
      - fvm::Sp(fvc::div(phi), U)
      + turbulence->divDevRhoReff(U)
    );

    UEqn().relax();

    eqnResidual = solve
    (
        UEqn() == -fvc::grad(pd) - fvc::grad(rho)*gh
    ).initialResidual();

    maxResidual = max(eqnResidual, maxResidual);

