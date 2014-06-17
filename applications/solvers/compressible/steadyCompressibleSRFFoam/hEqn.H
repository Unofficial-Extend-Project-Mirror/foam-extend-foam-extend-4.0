{
    // Solve the enthalpy equation
    T.storePrevIter();

    surfaceScalarField faceU = phi/fvc::interpolate(rho);

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
      - (turbulence->devRhoReff() && fvc::grad(Urel))
    );

    hEqn.relax();

    eqnResidual = hEqn.solve().initialResidual();
    maxResidual = max(eqnResidual, maxResidual);

    // Bound the enthalpy using TMin and TMax
    volScalarField Cp = thermo.Cp();

    h = Foam::min(h, TMax*Cp);
    h = Foam::max(h, TMin*Cp);
    h.correctBoundaryConditions();

    thermo.correct();
}
