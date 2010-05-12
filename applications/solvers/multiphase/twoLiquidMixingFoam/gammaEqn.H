{
    fvScalarMatrix gammaEqn
    (
        fvm::ddt(gamma)
      + fvm::div(phi, gamma)
      - fvm::laplacian(Dab, gamma)
    );

    gammaEqn.solve();

    rhoPhi = gammaEqn.flux()*(rho1 - rho2) + phi*rho2;
    rho = gamma*rho1 + (scalar(1) - gamma)*rho2;

    Info<< "Phase 1 volume fraction = "
        << gamma.weightedAverage(mesh.V()).value()
        << "  Min(gamma) = " << min(gamma).value()
        << "  Max(gamma) = " << max(gamma).value()
        << endl;
}
