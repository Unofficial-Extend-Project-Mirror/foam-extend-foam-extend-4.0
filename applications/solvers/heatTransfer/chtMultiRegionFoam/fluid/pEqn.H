{
    bool closedVolume = false;

    rhof[i] = thermof[i].rho();

    volScalarField rUA = 1.0/UEqn().A();
    Uf[i] = rUA*UEqn().H();

    phif[i] =
        fvc::interpolate(rhof[i])
       *(
            (fvc::interpolate(Uf[i]) & fluidRegions[i].Sf())
          + fvc::ddtPhiCorr(rUA, rhof[i], Uf[i], phif[i])
        )
      - fvc::interpolate(rhof[i]*rUA*ghf[i])
       *fvc::snGrad(rhof[i])
       *fluidRegions[i].magSf();

    // Solve pressure difference
#   include "pdEqn.H"

    // Solve continuity
#   include "rhoEqn.H"

    // Update pressure field (including bc)
    thermof[i].p() == pdf[i] + rhof[i]*ghf[i] + pRef;
    DpDtf[i] = fvc::DDt
    (
        surfaceScalarField("phiU", phif[i]/fvc::interpolate(rhof[i])),
        thermof[i].p()
    );

    // Update continuity errors
    compressibleContinuityErrors(cumulativeContErr, rhof[i], thermof[i]);

    // Correct velocity field
    Uf[i] -= rUA*(fvc::grad(pdf[i]) + fvc::grad(rhof[i])*ghf[i]);
    Uf[i].correctBoundaryConditions();

    // For closed-volume cases adjust the pressure and density levels
    // to obey overall mass continuity
    if (closedVolume)
    {
        thermof[i].p() +=
            (
                dimensionedScalar
                (
                    "massIni",
                    dimMass,
                    initialMassf[i]
                )
              - fvc::domainIntegrate(thermof[i].psi()*thermof[i].p())
            )/fvc::domainIntegrate(thermof[i].psi());
        pdf[i] == thermof[i].p() - (rhof[i]*ghf[i] + pRef);
        rhof[i] = thermof[i].rho();
    }

    // Update thermal conductivity
    Kf[i] = thermof[i].Cp()*turb[i].alphaEff();
}
