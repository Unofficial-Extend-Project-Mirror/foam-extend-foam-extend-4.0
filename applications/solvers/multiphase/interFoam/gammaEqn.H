{
    word gammaScheme("div(phi,gamma)");
    word gammarScheme("div(phirb,gamma)");

    surfaceScalarField phic = mag(phi/mesh.magSf());
    phic = min(interface.cGamma()*phic, max(phic));
    surfaceScalarField phir = phic*interface.nHatf();

    for (int gCorr=0; gCorr<nGammaCorr; gCorr++)
    {
        surfaceScalarField phiGamma = 
            fvc::flux
            (
                phi,
                gamma,
                gammaScheme
            )
          + fvc::flux
            (
                -fvc::flux(-phir, scalar(1) - gamma, gammarScheme),
                gamma,
                gammarScheme
            );

        MULES::explicitSolve(gamma, phi, phiGamma, 1, 0);

        rhoPhi = phiGamma*(rho1 - rho2) + phi*rho2;
    }

    Info<< "Liquid phase volume fraction = "
        << gamma.weightedAverage(mesh.V()).value()
        << "  Min(gamma) = " << min(gamma).value()
        << "  Max(gamma) = " << max(gamma).value()
        << endl;
}
