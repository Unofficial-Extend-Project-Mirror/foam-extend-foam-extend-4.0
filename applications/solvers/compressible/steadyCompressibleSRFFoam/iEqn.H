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
        // u & gradP term (steady-state formulation)
      + fvm::SuSp((fvc::div(faceU, p, "div(U,p)") - p*fvc::div(faceU))/i, i)
     ==
        // ddt(p) term removed: steady-state.  HJ, 27/Apr/2010
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
