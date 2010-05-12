{
    word gammaScheme("div(phi,gamma)");
    word gammarScheme("div(phirb,gamma)");

    surfaceScalarField phir("phir", phic*interface.nHatf());

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

        Pair<tmp<volScalarField> > vDotAlphal =
            twoPhaseProperties->vDotAlphal();
        const volScalarField& vDotcAlphal = vDotAlphal[0]();
        const volScalarField& vDotvAlphal = vDotAlphal[1]();

        volScalarField Sp
        (
            IOobject
            (
                "Sp",
                runTime.timeName(),
                mesh
            ),
            vDotvAlphal - vDotcAlphal
        );

        volScalarField Su
        (
            IOobject
            (
                "Su",
                runTime.timeName(),
                mesh
            ),
            // Divergence term is handled explicitly to be
            // consistent with the explicit transport solution
            divU*gamma
          + vDotcAlphal
        );

        //MULES::explicitSolve(gamma, phi, phiGamma, 1, 0);
        //MULES::explicitSolve(oneField(), gamma, phi, phiGamma, Sp, Su, 1, 0);
        MULES::implicitSolve(oneField(), gamma, phi, phiGamma, Sp, Su, 1, 0);

        rhoPhi +=
            (runTime.deltaT()/totalDeltaT)
           *(phiGamma*(rho1 - rho2) + phi*rho2);
    }

    Info<< "Liquid phase volume fraction = "
        << gamma.weightedAverage(mesh.V()).value()
        << "  Min(gamma) = " << min(gamma).value()
        << "  Max(gamma) = " << max(gamma).value()
        << endl;
}
