    // Solve the Momentum equation
    tmp<fvVectorMatrix> UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(phi, U)
      + turb.divDevRhoReff()
    );

    UEqn().relax();

    if (momentumPredictor)
    {
         solve
         (
            UEqn()
         ==
            fvc::reconstruct
            (
                fvc::interpolate(rho)*(g & mesh.Sf())
              - fvc::snGrad(p)*mesh.magSf()
            )
        );
    }
