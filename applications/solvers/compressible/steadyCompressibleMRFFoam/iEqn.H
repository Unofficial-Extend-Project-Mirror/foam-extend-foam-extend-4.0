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
        // u & gradP term (steady-state formulation)
      + fvm::SuSp((fvc::div(faceU, p, "div(U,p)") - fvc::div(faceU)*p)/i, i)
     ==
        // Viscous heating: note sign (devRhoReff has a minus in it)
      - (turbulence->devRhoReff() && fvc::grad(Urel))
    );

    iEqn.relax();

    eqnResidual = iEqn.solve().initialResidual();
    maxResidual = max(eqnResidual, maxResidual);

    // Calculate enthalpy out of rothalpy
    h = i + 0.5*magSqr(Urot);
    h.correctBoundaryConditions();

    // Bound the enthalpy using TMin and TMax
    volScalarField Cp = thermo.Cp();

    h = Foam::min(h, TMax*Cp);
    h = Foam::max(h, TMin*Cp);
    h.correctBoundaryConditions();

    // Re-initialise rothalpy based on limited enthalpy
    i = h - 0.5*magSqr(Urot);

    thermo.correct();
    psis = thermo.psi()/thermo.Cp()*thermo.Cv();
}
