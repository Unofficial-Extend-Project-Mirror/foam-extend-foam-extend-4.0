    const dictionary& stressControl =
        mesh.solutionDict().subDict("stressedFoam");

    int nCorr(readInt(stressControl.lookup("nCorrectors")));
    scalar convergenceTolerance(readScalar(stressControl.lookup("U")));
