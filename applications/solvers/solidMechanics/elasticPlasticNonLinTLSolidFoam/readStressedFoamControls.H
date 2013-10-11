const dictionary& stressControl =
  mesh.solutionDict().subDict("stressedFoam");

int nCorr(readInt(stressControl.lookup("nCorrectors")));
scalar convergenceTolerance(readScalar(stressControl.lookup("DU")));
Switch aitkenRelax(stressControl.lookup("aitkenRelaxation"));
//Switch thirdOrderCorrection(stressControl.lookup("thirdOrderCorrection"));
int infoFrequency(readInt(stressControl.lookup("infoFrequency")));
Switch largeStrainOverRelax(stressControl.lookup("largeStrainOverRelax"));
Switch nonLinearSemiImplicit(stressControl.lookup("nonLinearSemiImplicit"));
