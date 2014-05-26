    // Solve the momentum equation
    Urel.storePrevIter();

    fvVectorMatrix UrelEqn
    (
        fvm::ddt(rho, Urel)
      + fvm::div(phi, Urel)
      + turbulence->divDevRhoReff(Urel)
      + rho*SRF->Su()
    );

    UrelEqn.relax();

    eqnResidual = solve
    (
        UrelEqn == -fvc::grad(p)
    ).initialResidual();

    maxResidual = max(eqnResidual, maxResidual);
