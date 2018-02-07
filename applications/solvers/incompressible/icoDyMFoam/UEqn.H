    // Time-derivative matrix
    fvVectorMatrix ddtUEqn(fvm::ddt(U));

    // Convection-diffusion matrix
    fvVectorMatrix HUEqn
    (
        fvm::div(phi, U)
      - fvm::laplacian(nu, U)
    );

    if (piso.momentumPredictor())
    {
        solve(ddtUEqn + HUEqn == -fvc::grad(p));
    }
