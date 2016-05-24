    // Solve the momentum equation

    tmp<fvVectorMatrix> HUEqn
    (
        fvm::div(phi, U)
      + turbulence->divDevReff(U)
    );

    // Get under-relaxation factor and under-relax the diagonal
    const scalar UUrf = mesh.solutionDict().relaxationFactor(U.name());

    // Momentum solution
    eqnResidual = solve
    (
        relax(HUEqn(), UUrf)
     ==
        -fvc::grad(p)
    ).initialResidual();

    maxResidual = max(eqnResidual, maxResidual);
