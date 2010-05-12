    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(phi, U)
      + turbulence->divDevRhoReff(U)
     ==
        thermoCloud1.SU1()
      + kinematicCloud1.SU1()
      + rho.dimensionedInternalField()*g
    );

    UEqn.relax();

    if (momentumPredictor)
    {
        solve(UEqn == -fvc::grad(p));
    }
