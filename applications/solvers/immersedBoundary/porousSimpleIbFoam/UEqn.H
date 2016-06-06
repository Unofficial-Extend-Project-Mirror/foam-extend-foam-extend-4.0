    // Solve the Momentum equation
    tmp<fvVectorMatrix> UEqn
    (
        fvm::div(phi, U)
      + turbulence->divDevReff()
    );

    UEqn().relax();

    pZones.addResistance(UEqn());

    solve(UEqn() == -fvc::grad(p));
