{
    // Solve the enthalpy equation
    T.storePrevIter();

    // Calculate face velocity from flux
    surfaceScalarField faceU
    (
        "faceU",
        phi/fvc::interpolate(rho)
    );

    fvScalarMatrix hEqn
    (
        fvm::ddt(rho, h)
      + fvm::div(phi, h)
      - fvm::laplacian(turbulence->alphaEff(), h)
     ==
        fvc::div(faceU + mrfZones.fluxCorrection(), p, "div(U,p)")
        // Note: div flux correction is zero so there is no need to
        // carry it.  HJ, 4/Dec/2015
      - p*fvc::div(faceU)
        // Viscous heating: note sign (devRhoReff has a minus in it)
      - (turbulence->devRhoReff() && fvc::grad(U))
    );

    hEqn.relax();
    hEqn.solve();

    // Bounding of enthalpy taken out
    thermo.correct();
    psis = thermo.psi()/thermo.Cp()*thermo.Cv();
}
