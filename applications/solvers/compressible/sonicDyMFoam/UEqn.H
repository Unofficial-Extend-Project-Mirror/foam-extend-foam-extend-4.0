    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(phi, U)
      + turbulence->divDevRhoReff()
    );

    UEqn.relax
    (
        mesh.solutionDict().equationRelaxationFactor
        (
            U.select(pimple.finalIter())
        )
    );

    solve(UEqn == -fvc::grad(p));
