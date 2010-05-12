    surfaceScalarField gammaf = fvc::interpolate(gamma);

    surfaceScalarField muf =
        fvc::interpolate(rho)*(gammaf*nug + (1.0 - gammaf)*nul);

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(rhoPhi, U)
      - fvm::laplacian(muf, U)
      - (fvc::grad(U) & fvc::grad(muf))
    );
