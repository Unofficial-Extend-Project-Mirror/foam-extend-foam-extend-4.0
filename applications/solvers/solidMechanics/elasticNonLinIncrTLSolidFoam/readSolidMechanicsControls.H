    const dictionary& stressControl =
        mesh.solutionDict().subDict("solidMechanics");

    int nCorr(readInt(stressControl.lookup("nCorrectors")));
    scalar convergenceTolerance(readScalar(stressControl.lookup("DU")));
