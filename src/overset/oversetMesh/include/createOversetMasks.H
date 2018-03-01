//
// createOversetMask.H
// ~~~~~~~~~~~~
    const oversetMesh& om = oversetMesh::New(mesh);

    const volScalarField& cellOversetMask = om.gamma();

    const surfaceScalarField& faceOversetMask = om.sGamma();

    if (runTime.outputTime())
    {
#       include "writeOversetMasks.H"
    }
