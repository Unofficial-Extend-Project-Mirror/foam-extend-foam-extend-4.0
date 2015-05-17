const dictionary& stressControl =
  mesh.solutionDict().subDict("solidMechanics");

int nCorr(readInt(stressControl.lookup("nCorrectors")));
scalar convergenceTolerance(readScalar(stressControl.lookup("U")));
Switch predictor(stressControl.lookup("predictor"));
Switch dynamicTimeStep(stressControl.lookup("dynamicTimeStep"));
scalar deltaTmin(readScalar(stressControl.lookup("deltaTmin")));
scalar deltaTmax(readScalar(stressControl.lookup("deltaTmax")));
int infoFrequency(readInt(stressControl.lookup("infoFrequency")));
Switch aitkenRelax(stressControl.lookup("aitkenRelaxation"));
Switch relaxEqn(stressControl.lookup("relaxEqn"));

if(relaxEqn && solidInterfaceCorr)
{
    FatalError << "relaxEqn and solidInterface may not be used concurrently"
        << exit(FatalError);
}
