{
    // Create relative velocity
    Urel == U;
    mrfZones.relativeVelocity(Urel);

    // Bound the relative velocity to preserve the i to h conversion bound
    // HJ, 22/Mar/2017
    volScalarField magUrel = mag(Urel);

    if (max(magUrel) > UMax)
    {
        volScalarField UrelLimiter = pos(magUrel - UMax)*UMax/(magUrel + smallU)
            + neg(magUrel - UMax);
        UrelLimiter.max(scalar(0));
        UrelLimiter.min(scalar(1));

        Urel *= UrelLimiter;
        Urel.correctBoundaryConditions();
    }

    // Create rotational velocity (= omega x r)
    Urot == U - Urel;

    fvScalarMatrix iEqn
    (
        fvm::div(phi, i)
      + fvm::SuSp(-fvc::div(phi), i)
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
