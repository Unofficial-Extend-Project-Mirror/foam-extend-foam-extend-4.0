surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("0", dimensionSet(1, 0, -1, 0, 0), 0)
);

{
    label nGammaCorr
    (
        readLabel(piso.lookup("nGammaCorr"))
    );

    label nGammaSubCycles
    (
        readLabel(piso.lookup("nGammaSubCycles"))
    );

    surfaceScalarField phic = mag(phi/mesh.magSf());
    phic = min(interface.cGamma()*phic, max(phic));

    volScalarField divU = fvc::div(phi);

    dimensionedScalar totalDeltaT = runTime.deltaT();

    if (nGammaSubCycles > 1)
    {
        for
        (
            subCycle<volScalarField> gammaSubCycle(gamma, nGammaSubCycles);
            !(++gammaSubCycle).end();
        )
        {
#           include "gammaEqn.H"
        }
    }
    else
    {
#       include "gammaEqn.H"
    }

    if (nOuterCorr == 1)
    {
        interface.correct();
    }

    rho == gamma*rho1 + (scalar(1) - gamma)*rho2;
}
