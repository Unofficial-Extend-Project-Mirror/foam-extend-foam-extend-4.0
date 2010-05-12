    const dictionary& stressControl =
        stressMesh.solutionDict().subDict("stressedFoam");

    int nCorr(readInt(stressControl.lookup("nCorrectors")));
    scalar convergenceTolerance(readScalar(stressControl.lookup("U")));
