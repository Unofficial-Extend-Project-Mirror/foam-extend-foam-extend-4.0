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
        // ddt(p) term removed: steady-state.  HJ, 27/Apr/2010
        fvc::div(faceU, p, "div(U,p)")
      - p*fvc::div(faceU)
        // Viscous heating: note sign (devRhoReff has a minus in it)
      - (turbulence->devRhoReff() && fvc::grad(U))
    );

    hEqn.relax();

    eqnResidual = hEqn.solve().initialResidual();
    maxResidual = max(eqnResidual, maxResidual);

    // Bounding of enthalpy taken out
    thermo.correct();
}
