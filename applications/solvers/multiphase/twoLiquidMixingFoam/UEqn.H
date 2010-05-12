    surfaceScalarField muf = twoPhaseProperties.muf();

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(rhoPhi, U)
      - fvm::laplacian(muf, U)
        //- (fvc::grad(U) & fvc::grad(muf))
      - fvc::div(muf*(mesh.Sf() & fvc::interpolate(fvc::grad(U)().T())))
    );

    if (momentumPredictor)
    {
        solve
        (
            UEqn
         ==
           -fvc::reconstruct
            (
                mesh.magSf()*(fvc::snGrad(pd) + ghf*fvc::snGrad(rho))
            )
        );
    }
