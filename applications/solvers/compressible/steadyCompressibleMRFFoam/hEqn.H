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
        // Note: potential issue with reconstructed relative velocity.
        // HJ, 12/Dec/2009
//         fvc::div(faceU, p, "div(U,p)")
        fvc::div(faceU + mrfZones.fluxCorrection(), p, "div(U,p)")
      - p*fvc::div(faceU)
        // Viscous heating: note sign (devRhoReff has a minus in it)
      - (turbulence->devRhoReff() && fvc::grad(U))
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
