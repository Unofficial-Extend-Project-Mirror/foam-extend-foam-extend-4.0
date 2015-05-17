//- write force displacement to file
if(historyPatchID != -1)
{
    Info << "Writing disp-force to file for patch " << historyPatchName
        << endl;

    vector disp = gAverage(U.boundaryField()[historyPatchID]);

    //- for small strain or moving mesh
    vector force =
        gSum(mesh.boundary()[historyPatchID].Sf() & sigma.boundaryField()[historyPatchID]);

    //- for large strain total lagrangian
    // tensorField F = I + gradU.boundaryField()[historyPatchID];
    // vector force = gSum(mesh.Sf().boundaryField()[historyPatchID] & (sigma.boundaryField()[historyPatchID] & F));

    Info << "Writing strain-stress to file for patch " << historyPatchName
        << endl;

    // average stress strain
    symmTensor stress = gAverage(sigma.boundaryField()[historyPatchID]);
    symmTensor strain = gAverage(epsilon.boundaryField()[historyPatchID]);


    // write to file
    if(Pstream::master())
    {
        OFstream& forceDispFile = *forceFilePtr;
        label width = 20;
        forceDispFile << disp.x() << " " << disp.y() << " " << disp.z();
        forceDispFile.width(width);
        forceDispFile << force.x() << " " << force.y() << " " << force.z()
            << endl;

        OFstream& stressStrainFile = *stressFilePtr;
        stressStrainFile << strain.xx() << " " << strain.xy() << " " << strain.xz() << " "
            << strain.yy() << " " << strain.yz() << " " << strain.zz();
        stressStrainFile.width(width);
        stressStrainFile << stress.xx() << " " << stress.xy() << " " << stress.xz() << " "
            << stress.yy() << " " << stress.yz() << " " << stress.zz()
            << endl;
    }
}
