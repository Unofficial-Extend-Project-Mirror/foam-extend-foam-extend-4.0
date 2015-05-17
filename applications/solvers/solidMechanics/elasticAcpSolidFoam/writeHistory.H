//- write force displacement to file
if(historyPatchID != -1)
{
    Info << "Writing disp and force of patch "<<historyPatchName<<" to file"
        << endl;

    //- for small strain or moving mesh
    vector force = gSum(mesh.boundary()[historyPatchID].Sf() & sigma.boundaryField()[historyPatchID]);

    //- for large strain total lagrangian
    // tensorField F = I + gradU.boundaryField()[historyPatchID];
    // vectorField totalForce = mesh.Sf().boundaryField()[historyPatchID] & (sigma.boundaryField()[historyPatchID] & F);

    //vector force = sum( totalForce );

    //scalar force = force[vector::Y];

    //- patchIntegrate utility integrates it this way but this is worng because the sigma tensor should
    //- be dotted with the surface normal to give the actual traction/force
    //- you cannot just take the component of the sigma tensor
    //scalar forcePatchIntegrateMethod = gSum(
    //     mesh.magSf().boundaryField()[historyPatchID]*
    //     sigma.boundaryField()[historyPatchID].component(symmTensor::XY)
    //);

    vector avDisp = gAverage(U.boundaryField()[historyPatchID]);

    //- write to file
    if(Pstream::master())
    {
        OFstream& forceDispFile = *filePtr;
        forceDispFile << avDisp.x() << " " << avDisp.y() << " " << avDisp.z() << " "
            << force.x() << " " << force.y() << " " << force.z() << endl;
    }
}
