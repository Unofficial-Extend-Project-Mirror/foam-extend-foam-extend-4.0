{
    fvScalarMatrix eEqn
    (
        fvm::ddt(rho, e)
      + fvm::div(phi, e)
      - fvm::laplacian(turbulence->alphaEff(), e)
     ==
      - p*fvc::div(phi/fvc::interpolate(rho))
//       - fvm::SuSp
//         (
//             p*fvc::div(phi/fvc::interpolate(rho))/e,
//             e
//         )
        // viscous heating?
    );

    if (oCorr == nOuterCorr - 1)
    {
        if (mesh.solutionDict().relax("eFinal"))
        {
            eEqn.relax(mesh.solutionDict().relaxationFactor("eFinal"));
        }
        else
        {
            eEqn.relax(1);
        }
    }
    else
    {
        eEqn.relax();
    }

    eEqn.solve();

    thermo.correct();

    // Recalculate compressibility
    psis = thermo.psi()/thermo.Cp()*thermo.Cv();

    // Recalculate density
    rho = thermo.rho();
    rho.correctBoundaryConditions();
}
