{
    // Solve the enthalpy equation in total enthalpy formulation (see K)

    T.storePrevIter();

    K = 0.5*(magSqr(U));

    fvScalarMatrix hEqn
    (
        fvm::div(phi, h)
      + fvm::SuSp(-fvc::div(phi), h)
      + fvc::div(phi, K)
      - fvm::laplacian(turbulence->alphaEff(), h)
     ==
        // Viscous heating: note sign (devRhoReff has a minus in it)
      - (turbulence->devRhoReff() && fvc::grad(U))
    );

    hEqn.relax();
    hEqn.solve();

    // Bounding of enthalpy taken out
    thermo.correct();
    psis = thermo.psi()/thermo.Cp()*thermo.Cv();
}
