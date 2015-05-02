{
    fvScalarMatrix eEqn
    (
        fvm::ddt(rho, e)
      + fvm::div(phi, e)
      - fvm::laplacian(turbulence->alphaEff(), e)
     ==
      - fvm::SuSp(pDivU/e, e)
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

    // Bound the energy using TMin and TMax
    {
        dimensionedScalar Tstd("Tstd", dimTemperature, specie::Tstd());

        volScalarField Cv = thermo.Cv();
        volScalarField R = thermo.Cp() - Cv;

        e = Foam::min(e, TMax*Cv + R*Tstd);
        e = Foam::max(e, TMin*Cv + R*Tstd);
        e.correctBoundaryConditions();
    }

    thermo.correct();
}
