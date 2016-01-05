    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(phi, U)
      + turbulence->divDevRhoReff(U)
    );

    if (oCorr == nOuterCorr - 1)
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

    solve(UEqn == -fvc::grad(p));
