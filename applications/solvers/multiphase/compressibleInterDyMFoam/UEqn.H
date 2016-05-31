    surfaceScalarField muEff
    (
        "muEff",
        twoPhaseProperties.muf()
      + fvc::interpolate(rho*turbulence->nut())
    );

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(rhoPhi, U)
      - fvm::laplacian(muEff, U)
      - (fvc::grad(U) & fvc::grad(muEff))
      - fvc::div(muEff*(mesh.Sf() & fvc::interpolate(T(fvc::grad(U)))))
    );

    UEqn.relax();

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                fvc::interpolate(rho)*(g & mesh.Sf())
              + (
                    fvc::interpolate(interface.sigmaK())*fvc::snGrad(alpha1)
                  - fvc::snGrad(p)
                ) * mesh.magSf()
            )
        );
    }
