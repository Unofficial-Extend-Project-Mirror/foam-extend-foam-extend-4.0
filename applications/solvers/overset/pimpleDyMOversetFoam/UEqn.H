    fvVectorMatrix UEqn
    (
        fvm::ddt(U)
      + fvm::div(phi, U)
      + turbulence->divDevReff()
    );

    UEqn.relax
    (
        mesh.solutionDict().equationRelaxationFactor
        (
            U.select(pimple.finalIter())
        )
    );

    if (pimple.momentumPredictor())
    {
        solve(UEqn == -fvc::grad(p));
    }
