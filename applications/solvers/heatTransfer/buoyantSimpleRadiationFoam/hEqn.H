{
    fvScalarMatrix hEqn
    (
        fvm::div(phi, h)
      - fvm::Sp(fvc::div(phi), h)
      - fvm::laplacian(turbulence->alphaEff(), h)
     ==
        fvc::div(phi/fvc::interpolate(rho)*fvc::interpolate(p))
      - p*fvc::div(phi/fvc::interpolate(rho))
      + radiation->Sh(thermo)
    );

    hEqn.relax();

    hEqn.solve();

    thermo.correct();

    radiation->correct();
}
