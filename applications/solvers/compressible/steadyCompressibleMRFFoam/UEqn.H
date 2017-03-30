    fvVectorMatrix UEqn
    (
        fvm::div(phi, U)
      + fvm::SuSp(-fvc::div(phi), U)
      + turbulence->divDevRhoReff()
    );

    // MRF: add Coriolis force
    mrfZones.addCoriolis(rho, UEqn);

    UEqn.relax();

    solve(UEqn == -fvc::grad(p));
