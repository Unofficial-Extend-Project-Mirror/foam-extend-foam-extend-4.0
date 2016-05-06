    fvVectorMatrix UEqn
    (
        fvm::ddt(U)
      + fvm::div(phi, U)
      + turbulence->divDevReff(U)
    );

    UEqn.relax
    (
        mesh.solutionDict().relaxationFactor(U.select(pimple.finalIter()))
    );

    if (pimple.momentumPredictor())
    {
        solve(UEqn == -fvc::grad(p));
    }
