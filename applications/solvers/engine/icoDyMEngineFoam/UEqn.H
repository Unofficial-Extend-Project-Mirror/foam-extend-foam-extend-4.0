    fvVectorMatrix UEqn
    (
        fvm::ddt(U)
      + fvm::div(phi, U)
      - fvm::laplacian(nu, U)
    );

    if (momentumPredictor)
    {
        solve(UEqn == -fvc::grad(p));
    }
