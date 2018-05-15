    // Solve the momentum equation
    U.storePrevIter();

    fvVectorMatrix UEqn
    (
        fvm::div(phi, U)
      + turbulence->divDevRhoReff()
    );

    UEqn.relax();

    solve(UEqn == -fvc::grad(p));
