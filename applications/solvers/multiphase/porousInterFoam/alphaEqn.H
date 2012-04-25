{
    // Creates the porosity field for MULES
    volScalarField porosity
    (
        IOobject
        (
            "porosity",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("NULL", dimless, 1.0),
        "zeroGradient"
    );

    forAll (pZones, zoneI)
    {
        const label zoneId = pZones[zoneI].zoneId();

        const labelList& cells= mesh.cellZones()[zoneId];

        const scalar& zonePorosity = pZones[zoneI].porosity();

        forAll (cells, cellI)
        {
            porosity[cells[cellI]] = zonePorosity;
        }
    }

    // Ordinary MULES except for the argument list to explicitSolve
    word alphaScheme("div(phi,alpha)");
    word alpharScheme("div(phirb,alpha)");

    surfaceScalarField phic = mag(phi/mesh.magSf());
    phic = min(interface.cAlpha()*phic, max(phic));
    surfaceScalarField phir = phic*interface.nHatf();

    for (int aCorr=0; aCorr<nAlphaCorr; aCorr++)
    {
        surfaceScalarField phiAlpha =
            fvc::flux
            (
                phi,
                alpha1,
                alphaScheme
            )
          + fvc::flux
            (
                -fvc::flux(-phir, scalar(1) - alpha1, alpharScheme),
                alpha1,
                alpharScheme
            );

        MULES::explicitSolve
        (
            porosity,
            alpha1,
            phi,
            phiAlpha,
            zeroField(),
            zeroField(),
            1,
            0
        );

        rhoPhi = phiAlpha*(rho1 - rho2) + phi*rho2;
    }

    // The weightedAverage will change even in a closed box,
    // however, the volume of water should
    // remain constant - hence both values are written in the log
    Info<< "Liquid phase volume fraction = "
        << alpha1.weightedAverage(mesh.V()).value()
        << "  Volume of water = "
        << gSum( alpha1.internalField()*porosity.internalField()*mesh.V()
        << "  Min(alpha1) = " << min(alpha1).value()
        << "  Max(alpha1) = " << max(alpha1).value()
        << endl;
}
