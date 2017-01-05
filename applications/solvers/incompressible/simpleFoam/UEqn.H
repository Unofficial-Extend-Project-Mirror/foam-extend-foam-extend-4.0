    // Solve the momentum equation

    tmp<fvVectorMatrix> UEqn
    (
        fvm::div(phi, U)
      + turbulence->divDevReff()
    );

    // Momentum solution
    solve(relax(UEqn()) == -fvc::grad(p));
