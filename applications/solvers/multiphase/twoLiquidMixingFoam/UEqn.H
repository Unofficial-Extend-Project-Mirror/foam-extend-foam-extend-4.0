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
    //- fvc::div(muEff*(fvc::interpolate(dev(fvc::grad(U))) & mesh.Sf()))
    );

    if (oCorr == nOuterCorr-1)
    {
        if (mesh.solutionDict().relax("UFinal"))
        {
            UEqn.relax(mesh.solutionDict().relaxationFactor("UFinal"));
        }
        else
        {
            UEqn.relax(1);
        }
    }
    else
    {
        UEqn.relax();
    }

    if (momentumPredictor)
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (- ghf*fvc::snGrad(rho) - fvc::snGrad(p_rgh))*mesh.magSf()
            ),
            mesh.solutionDict().solver(oCorr == nOuterCorr-1 ? "UFinal" : "U")
        );
    }
