OFstream* forceFilePtr(NULL);
OFstream* stressFilePtr(NULL);
word historyPatchName(mesh.solutionDict().subDict("solidMechanics").lookup("historyPatch"));
label historyPatchID = mesh.boundaryMesh().findPatchID(historyPatchName);
if(historyPatchID == -1)
{
    Warning << "history patch " << historyPatchName
        << " not found. Force-displacement will not be written"
        << endl;
}
else if(Pstream::master())
{
    fileName historyDir = "history";
    mkDir(historyDir);

    {
        fileName forceFileName(historyDir/"forceDisp_"+historyPatchName+".dat");
        Info << "\nForce-displacement for patch " << historyPatchName
            << " will be written to " << forceFileName
            << endl;
        forceFilePtr = new OFstream(forceFileName);
        OFstream& forceDispFile = *forceFilePtr;
        forceDispFile << "#Disp(mm)\tForce(N)" << endl;
    }

    {
        fileName stressFileName(historyDir/"stressStrain_"+historyPatchName+".dat");
        Info << "\nStress(Engineering Small Stress)-strain(Engineering Small Strain) for patch "
            << historyPatchName
            << " will be written to " << stressFileName
            << endl;
        stressFilePtr = new OFstream(stressFileName);
        OFstream& stressStrainFile = *stressFilePtr;
        stressStrainFile << "#Strain(-)\tStress(Pa)" << endl;
    }
}
