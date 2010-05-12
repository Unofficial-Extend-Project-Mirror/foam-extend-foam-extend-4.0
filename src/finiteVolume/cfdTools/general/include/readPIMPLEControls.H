    dictionary pimple = mesh.solutionDict().subDict("PIMPLE");

    int nOuterCorr(readInt(pimple.lookup("nOuterCorrectors")));
    int nCorr(readInt(pimple.lookup("nCorrectors")));

    int nNonOrthCorr = 0;
    if (pimple.found("nNonOrthogonalCorrectors"))
    {
        nNonOrthCorr = readInt(pimple.lookup("nNonOrthogonalCorrectors"));
    }

    bool momentumPredictor = true;
    if (pimple.found("momentumPredictor"))
    {
        momentumPredictor = Switch(pimple.lookup("momentumPredictor"));
    }

    bool transonic = false;
    if (pimple.found("transonic"))
    {
        transonic = Switch(pimple.lookup("transonic"));
    }
