{
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
       - (turbulence->devRhoReff() && fvc::grad(U))
    );

    iEqn.relax();

    iEqn.solve();

    // From rothalpy, calculate enthalpy after solution of rothalpy equation
    h = i + 0.5*(magSqr(Urot) - magSqr(Urel));
    h.correctBoundaryConditions();
    // Update thermo for new h
    thermo.correct();
    psis = thermo.psi()/thermo.Cp()*thermo.Cv();
}

