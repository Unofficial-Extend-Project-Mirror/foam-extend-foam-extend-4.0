    // Solve the momentum equation

    tmp<fvVectorMatrix> HUrelEqn
    (
        fvm::div(phi, Urel)
      + turbulence->divDevReff()
      + SRF->Su()
    );

    // Get under-relaxation factor
    const scalar UUrf =
        mesh.solutionDict().equationRelaxationFactor(Urel.name());

    // Momentum solution
    solve(relax(HUrelEqn(), UUrf) == -fvc::grad(p));
