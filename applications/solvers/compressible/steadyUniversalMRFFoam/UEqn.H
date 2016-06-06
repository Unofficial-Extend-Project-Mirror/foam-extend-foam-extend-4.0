    // Solve the momentum equation
    U.storePrevIter();

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(phi, U)
      + turbulence->divDevRhoReff()
    );

    // MRF: add Coriolis force
    mrfZones.addCoriolis(rho, UEqn);

    UEqn.relax();

    eqnResidual = solve
    (
        UEqn == -fvc::grad(p)
    ).initialResidual();

    maxResidual = max(eqnResidual, maxResidual);
