    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(phi, U)
      + turbulence->divDevRhoReff()
     ==
        thermoCloud1.SU()
      + kinematicCloud1.SU()
      + rho.dimensionedInternalField()*g
    );

    UEqn.relax();

    if (pimple.momentumPredictor())
    {
        solve(UEqn == -fvc::grad(p));
    }
