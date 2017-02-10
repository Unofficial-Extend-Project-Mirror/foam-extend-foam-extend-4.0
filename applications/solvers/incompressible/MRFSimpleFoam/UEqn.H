    // Solve the momentum equation

    tmp<fvVectorMatrix> HUEqn
    (
        fvm::div(phi, U)
      + turbulence->divDevReff()
    );
    mrfZones.addCoriolis(HUEqn());

    // Get under-relaxation factor
    const scalar UUrf = mesh.solutionDict().equationRelaxationFactor(U.name());

    // Momentum solution
    solve
    (
        relax(HUEqn(), UUrf)
     ==
        -fvc::grad(p)
    );
