{
    // Solve the rothalpy equation
    T.storePrevIter();

    // Create relative velocity
    Urel == U;
    mrfZones.relativeVelocity(Urel);

    // Create rotational velocity (= omega x r)
    Urot = U - Urel;

    // Calculate face velocity from absolute flux
    surfaceScalarField rhof = fvc::interpolate(rho);

    surfaceScalarField phiAbs
    (
        "phiAbs",
        phi
    );
    mrfZones.absoluteFlux(rhof, phiAbs);

    surfaceScalarField faceU("faceU", phiAbs/rhof);

    fvScalarMatrix iEqn
    (
        fvm::ddt(rho, i)
      + fvm::div(phi, i)
      - fvm::laplacian(turbulence->alphaEff(), i)
     ==
        // ddt(p) term removed: steady-state.  HJ, 27/Apr/2010
        fvc::div(faceU, p, "div(U,p)")
      - p*fvc::div(faceU)
        // Viscous heating: note sign (devRhoReff has a minus in it)
      - (turbulence->devRhoReff() && fvc::grad(Urel))
    );

    iEqn.relax();

    eqnResidual = iEqn.solve().initialResidual();
    maxResidual = max(eqnResidual, maxResidual);

    // Calculate enthalpy out of rothalpy
    h = i + 0.5*magSqr(Urot);
    h.correctBoundaryConditions();

    thermo.correct();
    psis = thermo.psi()/thermo.Cp()*thermo.Cv();
}
