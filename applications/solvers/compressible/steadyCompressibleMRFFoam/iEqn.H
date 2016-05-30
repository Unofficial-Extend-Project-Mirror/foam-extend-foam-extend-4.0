{
    // Solve the rothalpy equation

    // Create relative velocity
    Urel == U;
    mrfZones.relativeVelocity(Urel);

    // Create rotational velocity (= omega x r)
    Urot == U - Urel;

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

    iEqn.solve();

    // Calculate enthalpy out of rothalpy
    h = i + 0.5*(magSqr(Urot) - magSqr(Urel));
    h.correctBoundaryConditions();

    thermo.correct();
    psis = thermo.psi()/thermo.Cp()*thermo.Cv();
}
