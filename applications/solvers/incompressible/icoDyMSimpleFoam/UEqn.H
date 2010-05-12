    fvVectorMatrix UEqn
    (
        fvm::ddt(U)
      + fvm::div(phi, U)
      - fvm::laplacian(nu, U)
    );

    UEqn.relax();

    // Solve the momentum equation
    solve(UEqn == -fvc::grad(p));
