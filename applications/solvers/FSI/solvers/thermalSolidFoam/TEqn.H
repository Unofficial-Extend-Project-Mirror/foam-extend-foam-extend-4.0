{
    int iCorr = 0;
    scalar relResidual = GREAT;

    do
    {
        T.storePrevIter();

        fvScalarMatrix TEqn
        (
            rhoC*fvm::ddt(T)
         == fvm::laplacian(k, T, "laplacian(k,T)")
        );

        TEqn.solve();

        // Calc rel. residual
        {
            scalar maxDT =
                gMax(mag(T.internalField() - T.oldTime().internalField()));

            relResidual =
                gMax
                (
                    mag(T.internalField() - T.prevIter().internalField())
                   /(maxDT + SMALL)
                );
        }
        Info << "Relative temperature residual: " << relResidual << endl;
    }
    while((relResidual > 1e-6) && (++iCorr < 10));
}

volScalarField DT("DT", T - T0);
surfaceScalarField DTf("DTf", fvc::interpolate(DT));
