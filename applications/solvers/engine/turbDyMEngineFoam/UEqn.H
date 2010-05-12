    fvVectorMatrix UEqn
    (
        fvm::ddt(U)
      + fvm::div(phi, U)
      + turbulence->divDevReff(U)
    );

    if (momentumPredictor)
    {
        solve(UEqn == -fvc::grad(p));
    }
