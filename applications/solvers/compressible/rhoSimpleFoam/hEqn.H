{
    fvScalarMatrix hEqn
    (
        fvm::div(phi, h)
      - fvm::Sp(fvc::div(phi), h)
      - fvm::laplacian(turbulence->alphaEff(), h)
     ==
        fvc::div(phi/fvc::interpolate(rho), rho/psi, "div(U,p)")
      - (rho/psi)*fvc::div(phi/fvc::interpolate(rho))
    );

    hEqn.relax();

    hEqn.solve();

    thermo.correct();
}
