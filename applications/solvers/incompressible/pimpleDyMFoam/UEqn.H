    // Convection-diffusion matrix
    fvVectorMatrix HUEqn
    (
        fvm::div(phi, U)
      + turbulence->divDevReff(U)
    );

    // Time derivative matrix
    fvVectorMatrix ddtUEqn(fvm::ddt(U));

    // Get under-relaxation factor
    scalar UUrf = mesh.solutionDict().relaxationFactor(U.name());

    if (oCorr == nOuterCorr - 1)
    {
        if (mesh.solutionDict().relax("UFinal"))
        {
            UUrf = mesh.solutionDict().relaxationFactor("UFinal");
        }
        else
        {
            UUrf = 1;
        }
    }

    // Solve momentum predictor
    solve
    (
        ddtUEqn
      + relax(HUEqn, UUrf)
     ==
      - fvc::grad(p)
    );
