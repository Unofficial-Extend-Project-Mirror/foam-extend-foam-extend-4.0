    // Solve the Momentum equation

    tmp<fvVectorMatrix> UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(phi, U)
      + turbulence->divDevRhoReff(U)
    );

//if (oCorr != nOuterCorr-1)
    {
        UEqn().relax();
    }

    if (momentumPredictor)
    {
        if (oCorr == nOuterCorr-1)
        {
            solve(UEqn() == -fvc::grad(p), mesh.solver("UFinal"));
        }
        else
        {
            solve(UEqn() == -fvc::grad(p));
        }
    }
