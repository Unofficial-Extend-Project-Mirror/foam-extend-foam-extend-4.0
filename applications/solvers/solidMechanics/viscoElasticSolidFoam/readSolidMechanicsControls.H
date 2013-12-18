    const dictionary& stressControl =
        mesh.solutionDict().subDict("solidMechanics");

    int nCorr(readInt(stressControl.lookup("nCorrectors")));
    scalar convergenceTolerance(readScalar(stressControl.lookup("DU")));

//     componentReferenceList cr
//     (
//         stressControl.lookup("componentReference"),
//         componentReference::iNew(mesh)
//     );


//     dimensionedScalar K
//     (
//         mesh.solutionDict().subDict("solidMechanics").lookup("K")
//     );

int infoFrequency(readInt(stressControl.lookup("infoFrequency")));

dimensionedVector g("gravity", dimLength/(dimTime*dimTime), stressControl.lookup("gravity"));
