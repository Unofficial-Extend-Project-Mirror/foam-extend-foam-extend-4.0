    // Solve the Momentum equation

    tmp<fvVectorMatrix> UEqn
    (
        fvm::div(phi, U)
      + turbulence->divDevRhoReff()
    );

    UEqn().relax();

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
