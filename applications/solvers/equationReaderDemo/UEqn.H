    // Solve the Momentum equation

    tmp<fvVectorMatrix> UEqn
    (
        fvm::div(phi, U)
      + turbulence->divDevReff(U)
    );

    UEqn().relax();

    solve(UEqn() == -fvc::grad(p));
