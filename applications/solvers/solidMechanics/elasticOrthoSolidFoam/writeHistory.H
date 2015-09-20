//- write force displacement to file
if(historyPatchID != -1)
{
    Info << "Writing disp and force of patch "<<historyPatchName<<" to file"
        << endl;

    //- for small strain or moving mesh
    vector force = gSum(mesh.boundary()[historyPatchID].Sf() & sigma.boundaryField()[historyPatchID]);

    vector avDisp = gAverage(U.boundaryField()[historyPatchID]);

    //- write to file
    if(Pstream::master())
    {
        OFstream& forceDispFile = *filePtr;
        forceDispFile << avDisp.x() << " " << avDisp.y() << " " << avDisp.z() << " "
            << force.x() << " " << force.y() << " " << force.z() << endl;
    }
}
