{
    word alphaScheme("div(phi,alpha)");
    word alpharScheme("div(phirb,alpha)");

    // New formulation of phir: Co condition
    surfaceScalarField phir
    (
        "phir",
        faceIbMask*interface.cAlpha()*interface.nHatf()*
        min
        (
            scalar(0.5)/    // This is ref compression Co number for cAlpha = 1
            (
                mesh.surfaceInterpolation::deltaCoeffs()*
                mesh.time().deltaT()
            ),
            mag(phi/mesh.magSf())
          + dimensionedScalar("small", dimVelocity, 1e-3)
        )
    );

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

    alpha1Eqn.relax();
    alpha1Eqn.solve();

    Info<< "Liquid phase volume fraction = "
        << alpha1.weightedAverage(mesh.V()).value()
        << "  Min(alpha1) = " << min(cellIbMask*alpha1).value()
        << "  Max(alpha1) = " << max(cellIbMask*alpha1).value()
        << endl;

    rhoPhi = faceIbMask*(alpha1Eqn.flux()*(rho1 - rho2) + phi*rho2);

    // Limit alpha to handle IB cells
//     alpha1.max(0);
//     alpha1.min(1);

    // Update of interface and density
    interface.correct();
    twoPhaseProperties.correct();

    // Calculate density using limited alpha1
    rho = twoPhaseProperties.rho();
    rho.correctBoundaryConditions();
}
