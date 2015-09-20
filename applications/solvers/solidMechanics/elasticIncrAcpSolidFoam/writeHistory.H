//- write force displacement to file
if(historyPatchID != -1)
{
    Info << "Found patch "<<historyPatchName<<", writing y force and displacement to file"
        << endl;

    //- calculate force in specified direction on topClamp patch
    vector direction(0, 1, 0);

    //- for small strain or moving mesh
    scalar force = gSum
    (
        direction &
        (mesh.boundary()[historyPatchID].Sf() & sigma.boundaryField()[historyPatchID])
    );

    //- for large strain total lagrangian
    // tensorField F = I + gradU.boundaryField()[historyPatchID];
    // vectorField totalForce = mesh.Sf().boundaryField()[historyPatchID] & (sigma.boundaryField()[historyPatchID] & F);

    //vector force = sum( totalForce );

    //scalar force = force[vector::Y];

    //- patchIntegrate utility integrates it this way but this is worng because the sigma tensor should
    //- be dotted with the surface normal to give the actual traction/force
    //- you cannot just take the component of the sigma tensor
    //scalar forcePatchIntegrateMethod = gSum(
    //                        mesh.magSf().boundaryField()[historyPatchID]
    //                        *sigma.boundaryField()[historyPatchID].component(symmTensor::XY)
    //                        );

    scalar disp = max(U.boundaryField()[historyPatchID].component(vector::Y));

    //- write to file
    if(Pstream::master())
    {
        OFstream& forceDispFile = *filePtr;
        forceDispFile << disp << "\t" << force << endl;
    }
}
