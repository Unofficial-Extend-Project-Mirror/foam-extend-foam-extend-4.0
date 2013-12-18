const dictionary& stressControl =
  mesh.solutionDict().subDict("solidMechanics");

int nCorr(readInt(stressControl.lookup("nCorrectors")));
scalar convergenceTolerance(readScalar(stressControl.lookup("DU")));

Switch aitkenRelax(stressControl.lookup("aitkenRelaxation"));

int infoFrequency(readInt(stressControl.lookup("infoFrequency")));
Switch nonLinearSemiImplicit(stressControl.lookup("nonLinearSemiImplicit"));
