const dictionary& stressControl =
  mesh.solutionDict().subDict("solidMechanics");

int nCorr(readInt(stressControl.lookup("nCorrectors")));

scalar convergenceToleranceU(readScalar(stressControl.lookup("U")));

scalar convergenceToleranceT(readScalar(stressControl.lookup("U")));

int infoFrequency(readInt(stressControl.lookup("infoFrequency")));

Switch aitkenRelax(stressControl.lookup("aitkenRelaxation"));
