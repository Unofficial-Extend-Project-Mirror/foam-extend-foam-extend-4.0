const dictionary& stressControl =
  mesh.solutionDict().subDict("solidMechanics");

int nCorr(readInt(stressControl.lookup("nCorrectors")));
scalar convergenceTolerance(readScalar(stressControl.lookup("DU")));

// componentReferenceList cr
// (
//  stressControl.lookup("componentReference"),
//  componentReference::iNew(mesh)
//  );

Switch dynamicTimeStep(stressControl.lookup("dynamicTimeStep"));
scalar deltaTmin(readScalar(stressControl.lookup("deltaTmin")));
scalar deltaTmax(readScalar(stressControl.lookup("deltaTmax")));

int infoFrequency(readInt(stressControl.lookup("infoFrequency")));
//Info << "\nFrequency of info statements is " << infoFrequency << endl;

Switch aitkenRelax(stressControl.lookup("aitkenRelaxation"));
Switch predictor(stressControl.lookup("predictor"));
