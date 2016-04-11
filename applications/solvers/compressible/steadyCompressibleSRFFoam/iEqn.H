{
    // Solve the rothalpy equation
    T.storePrevIter();

    // Calculate face velocity from flux
    surfaceScalarField faceU
    (
        "faceU",
        phi/fvc::interpolate(rho) + (SRF->faceU() & mesh.Sf())
    );

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
    volVectorField Urot("Urot", SRF->U());

    h = i + 0.5*magSqr(Urot);
    h.correctBoundaryConditions();

    // Re-initialise rothalpy based on limited enthalpy
    i = h - 0.5*magSqr(Urot);

    // Bounding of enthalpy taken out
    thermo.correct();
    psis = thermo.psi()/thermo.Cp()*thermo.Cv();
}
