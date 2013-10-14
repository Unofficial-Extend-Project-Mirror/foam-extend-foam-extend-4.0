const dictionary& stressControl =
  mesh.solutionDict().subDict("solidMechanics");

int nCorr(readInt(stressControl.lookup("nCorrectors")));
scalar convergenceTolerance(readScalar(stressControl.lookup("U")));
Switch predictor(stressControl.lookup("predictor"));
int infoFrequency(readInt(stressControl.lookup("infoFrequency")));
Switch aitkenRelax(stressControl.lookup("aitkenRelaxation"));
//Switch thirdOrderCorrection(stressControl.lookup("thirdOrderCorrection"));
//Switch relaxEqn(stressControl.lookup("relaxEquation"));

// if(relaxEqn && solidInterfaceCorr)
//   {
//     FatalError << "relaxEqn and solidInterface may not be used concurrently"
//                << exit(FatalError);
//   }
