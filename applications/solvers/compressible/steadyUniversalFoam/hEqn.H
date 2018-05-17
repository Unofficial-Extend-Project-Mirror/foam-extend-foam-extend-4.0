{
    // Solve the enthalpy equation

    // Calculate face velocity from flux
    surfaceScalarField faceU
    (
        "faceU",
        phi/fvc::interpolate(rho)
    );

    fvScalarMatrix hEqn
    (
        fvm::div(phi, h)
      + fvm::SuSp(-fvc::div(phi), h)
      - fvm::laplacian(turbulence->alphaEff(), h)
     ==
        fvc::div(faceU, p, "div(U,p)")
      - p*fvc::div(faceU)
        // ddt(p) term removed: steady-state.  HJ, 27/Apr/2010
        // Viscous heating: note sign (devRhoReff has a minus in it)
      - (turbulence->devRhoReff() && fvc::grad(U))
    );

    hEqn.relax();
    hEqn.solve();

    // Bounding of enthalpy taken out
    thermo.correct();

    if (compressible)
    {
        psis = thermo.psi()/thermo.Cp()*thermo.Cv();
    }
}
