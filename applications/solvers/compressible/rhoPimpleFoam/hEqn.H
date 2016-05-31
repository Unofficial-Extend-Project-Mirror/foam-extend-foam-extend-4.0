{
    fvScalarMatrix hEqn
    (
        fvm::ddt(rho, h)
      + fvm::div(phi, h)
      - fvm::laplacian(turbulence->alphaEff(), h)
     ==
        DpDt
    );

    hEqn.relax
    (
        mesh.solutionDict().equationRelaxationFactor(h.select(pimple.finalIter()))
    );
    hEqn.solve(mesh.solutionDict().solver((h.select(pimple.finalIter()))));

    thermo.correct();
}
