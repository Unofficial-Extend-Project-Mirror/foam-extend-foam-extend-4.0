{
    // Solve the rothalpy equation
    T.storePrevIter();

    // Create relative velocity
    Urel == U;
    mrfZones.relativeVelocity(Urel);

    // Create rotational velocity (= omega x r)
    Urot = U - Urel;

    fvScalarMatrix iEqn
    (
        fvm::ddt(rho, i)
      + fvm::div(phi, i)
      - fvm::laplacian(turbulence->alphaEff(), i)
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
}
