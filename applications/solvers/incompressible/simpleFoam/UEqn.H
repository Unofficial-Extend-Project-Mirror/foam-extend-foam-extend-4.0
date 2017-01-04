    // Solve the momentum equation

    tmp<fvVectorMatrix> UEqn
    (
        fvm::div(phi, U)
      + turbulence->divDevReff()
    );

    // Relax the equation
    UEqn().relax();

    // Momentum solution
    solve(UEqn() == -fvc::grad(p));
