    // Solve the momentum equation
    Urel.storePrevIter();

    fvVectorMatrix UrelEqn
    (
        fvm::ddt(rho, Urel)
      + fvm::div(phi, Urel)
      + turbulence->divDevRhoReff()
      + rho*SRF->Su()
    );

    UrelEqn.relax();

    solve(UrelEqn == -fvc::grad(p));
