    // Momentum equation
    fvVectorMatrix UEqn
    (
        fvm::ddt(U)
      + fvm::div(phi, U)
      + turbulence->divDevReff(U)
    );

    // Add MRF and porous sources
    mrfZones.addCoriolis(UEqn);
    pZones.addResistance(UEqn);

    UEqn.relax();

    UpEqn.insertEquation(0, UEqn);
