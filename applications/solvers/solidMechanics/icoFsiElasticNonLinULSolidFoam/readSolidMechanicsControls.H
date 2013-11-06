    const dictionary& stressControl =
        stressMesh.solutionDict().subDict("solidMechanics");

    int nCorr(readInt(stressControl.lookup("nCorrectors")));
    scalar convergenceTolerance(readScalar(stressControl.lookup("DU")));
