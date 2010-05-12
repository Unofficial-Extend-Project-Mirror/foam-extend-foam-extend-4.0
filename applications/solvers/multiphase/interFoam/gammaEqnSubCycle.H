label nGammaCorr
(
    readLabel(piso.lookup("nGammaCorr"))
);

label nGammaSubCycles
(
    readLabel(piso.lookup("nGammaSubCycles"))
);

if (nGammaSubCycles > 1)
{
    dimensionedScalar totalDeltaT = runTime.deltaT();
    surfaceScalarField rhoPhiSum = 0.0*rhoPhi;

    for
    (
        subCycle<volScalarField> gammaSubCycle(gamma, nGammaSubCycles);
        !(++gammaSubCycle).end();
    )
    {
#       include "gammaEqn.H"
        rhoPhiSum += (runTime.deltaT()/totalDeltaT)*rhoPhi;
    }

    rhoPhi = rhoPhiSum;
}
else
{
#       include "gammaEqn.H"
}

interface.correct();

rho == gamma*rho1 + (scalar(1) - gamma)*rho2;
