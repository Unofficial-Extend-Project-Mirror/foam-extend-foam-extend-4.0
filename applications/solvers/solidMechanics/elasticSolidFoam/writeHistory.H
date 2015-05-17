//- write force displacement to file
if(historyPatchID != -1)
{
    Info << "Found patch "<<historyPatchName<<", writing y force and displacement to file"
        << endl;

    //- for small strain or moving mesh
    vector force =
        gSum(mesh.boundary()[historyPatchID].Sf() & sigma.boundaryField()[historyPatchID]);

    vector disp = gAverage(U.boundaryField()[historyPatchID]);

    //- write to file
    if(Pstream::master())
    {
        OFstream& forceDispFile = *filePtr;
        forceDispFile << disp.x() << " "
            << disp.y() << " "
            << disp.z() << " "
            << force.x() << " "
            << force.y() << " "
            << force.z() << endl;
    }
}
