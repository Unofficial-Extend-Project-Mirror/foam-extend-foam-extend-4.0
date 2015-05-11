{
    word alphaScheme("div(phi,alpha)");
    word alpharScheme("div(phirb,alpha)");

    surfaceScalarField phic = mag(phi/mesh.magSf());
    phic = min(interface.cAlpha()*phic, max(phic));
    surfaceScalarField phir = faceIbMask*phic*interface.nHatf();

    fvScalarMatrix alpha1Eqn
    (
        fvm::ddt(alpha1)
      + fvm::div(phi, alpha1, alphaScheme)
      + fvm::div
        (
            -fvc::flux(-phir, scalar(1) - alpha1, alpharScheme),
            alpha1,
            alpharScheme
        )
    );

    alpha1Eqn.solve();

    Info<< "Liquid phase volume fraction = "
        << alpha1.weightedAverage(mesh.V()).value()
        << "  Min(alpha1) = " << min(cellIbMask*alpha1).value()
        << "  Max(alpha1) = " << max(cellIbMask*alpha1).value()
        << endl;

    // Limit alpha to handle IB cells
    alpha1.max(0);
    alpha1.min(1);

    // Update of interface and density
    interface.correct();

    rho == alpha1*rho1 + (scalar(1) - alpha1)*rho2;
}
