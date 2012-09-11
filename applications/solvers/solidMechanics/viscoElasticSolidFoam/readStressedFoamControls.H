    const dictionary& stressControl =
        mesh.solutionDict().subDict("stressedFoam");

    int nCorr(readInt(stressControl.lookup("nCorrectors")));
    scalar convergenceTolerance(readScalar(stressControl.lookup("DU")));

//     componentReferenceList cr
//     (
//         stressControl.lookup("componentReference"),
//         componentReference::iNew(mesh)
//     );


//     dimensionedScalar K
//     (
//         mesh.solutionDict().subDict("stressedFoam").lookup("K")
//     );
