{
    // Solve the rothalpy equation

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

    T.storePrevIter();

    fvScalarMatrix iEqn
    (
        fvm::div(phi, i)
      + fvm::SuSp(-fvc::div(phi), i)
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
