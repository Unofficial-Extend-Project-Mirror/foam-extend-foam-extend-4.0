OFstream * filePtr(NULL);
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
    Info << "Force-displacement for patch " << historyPatchName
        << " will be written to forceDisp.dat"
        << endl;
    filePtr = new OFstream("forceDisp.dat");
    OFstream& forceDispFile = *filePtr;
    forceDispFile << "#Disp(mm)\tForce(N)" << endl;
}
