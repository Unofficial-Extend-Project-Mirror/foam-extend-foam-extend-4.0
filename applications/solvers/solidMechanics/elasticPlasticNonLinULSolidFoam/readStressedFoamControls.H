    const dictionary& stressControl =
        mesh.solutionDict().subDict("stressedFoam");

    int nCorr(readInt(stressControl.lookup("nCorrectors")));
    scalar convergenceTolerance(readScalar(stressControl.lookup("DU")));

//     scalar interfaceRelax(readScalar(stressControl.lookup("interfaceRelaxation")));

//     componentReferenceList cr
//     (
//         stressControl.lookup("componentReference"),
//         componentReference::iNew(mesh)
//    );
