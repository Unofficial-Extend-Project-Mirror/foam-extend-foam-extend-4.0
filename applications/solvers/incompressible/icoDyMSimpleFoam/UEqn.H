    // Time-derivative matrix
    fvVectorMatrix ddtUEqn(fvm::ddt(U));

    // Convection-diffusion matrix
    fvVectorMatrix HUEqn
    (
        fvm::div(phi, U)
      - fvm::laplacian(nu, U)
    );

    // Solve the relaxed momentum equation
    solve(relax(ddtUEqn + HUEqn) == -fvc::grad(p));
