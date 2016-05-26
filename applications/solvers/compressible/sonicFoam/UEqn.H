    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(phi, U)
      + turbulence->divDevRhoReff(U)
    );

    UEqn.relax
    (
        mesh.solutionDict().relaxationFactor(U.select(pimple.finalIter()))
    );

    solve(UEqn == -fvc::grad(p));
