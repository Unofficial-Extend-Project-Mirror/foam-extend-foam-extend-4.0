    const dictionary& simple = mesh.solutionDict().subDict("SIMPLE");

    int nTimeCorr = 0;

    if (simple.found("nTimeCorrectors"))
    {
        nTimeCorr = readInt(simple.lookup("nTimeCorrectors"));
    }

    int nNonOrthCorr = 0;

    if (simple.found("nNonOrthogonalCorrectors"))
    {
        nNonOrthCorr = readInt(simple.lookup("nNonOrthogonalCorrectors"));
    }
