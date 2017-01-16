    // Solve the momentum equation

    tmp<fvVectorMatrix> UrelEqn
    (
        fvm::div(phi, Urel)
      + turbulence->divDevReff()
      + SRF->Su()
    );

    // Momentum solution
    solve(relax(UrelEqn()) == -fvc::grad(p));
