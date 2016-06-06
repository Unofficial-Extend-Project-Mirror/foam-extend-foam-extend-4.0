    // Solve the momentum equation
    U.storePrevIter();

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(phi, U)
      + turbulence->divDevRhoReff()
    );

    UEqn.relax();

    eqnResidual = solve
    (
        UEqn == -fvc::grad(p)
    ).initialResidual();

    maxResidual = max(eqnResidual, maxResidual);
