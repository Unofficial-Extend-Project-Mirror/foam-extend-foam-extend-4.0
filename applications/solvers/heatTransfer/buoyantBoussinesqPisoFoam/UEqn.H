    // Solve the momentum equation

    fvVectorMatrix UEqn
    (
        fvm::ddt(U)
      + fvm::div(phi, U)
      + turbulence->divDevReff()
    );

    UEqn.relax();

    if (piso.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                    fvc::interpolate(rhok)*(g & mesh.Sf())
                  - fvc::snGrad(p)*mesh.magSf()
                )
            )
        );
    }
