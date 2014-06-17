// Loading patch
{
    label patchID = mesh.boundaryMesh().findPatchID("loading");

    if(patchID == -1)
    {
        FatalErrorIn(args.executable())
            << "Can't find patch with name loading" << abort(FatalError);
    }

    if
    (
        DU.boundaryField()[patchID].type()
     == directionMixedFvPatchVectorField::typeName
    )
    {
        directionMixedFvPatchVectorField& loadingPatchU =
            refCast<directionMixedFvPatchVectorField>
            (
                DU.boundaryField()[patchID]
            );

        vectorField n = mesh.boundary()[patchID].nf();

        vectorField traction = n*(n&(n&DSigma));

        scalarField lpMu = mu.boundaryField()[patchID];
        scalarField lpLambda = lambda.boundaryField()[patchID];

        tensorField lpGradDU = gradDU.boundaryField()[patchID];

        loadingPatchU.refGrad() =
        (
            traction
          - (n & (lpMu*lpGradDU.T() - (lpMu + lpLambda)*lpGradDU))
          - n*lpLambda*tr(lpGradDU)
        )
       /(2.0*lpMu + lpLambda);
    }
}

// Bottom patch
{
    label patchID = mesh.boundaryMesh().findPatchID("bottom");

    if(patchID == -1)
    {
        FatalErrorIn(args.executable())
            << "Can't find patch with name bottom" << abort(FatalError);
    }

    if
    (
        DU.boundaryField()[patchID].type()
     == directionMixedFvPatchVectorField::typeName
    )
    {
        directionMixedFvPatchVectorField& bottomPatchU =
            refCast<directionMixedFvPatchVectorField>
            (
                DU.boundaryField()[patchID]
            );

        vectorField n = mesh.boundary()[patchID].nf();

        vectorField traction = n*(n&(n&DSigma));

        scalarField bpMu = mu.boundaryField()[patchID];
        scalarField bpLambda = lambda.boundaryField()[patchID];

        tensorField bpGradDU = gradDU.boundaryField()[patchID];

        bottomPatchU.refGrad() =
        (
            traction
          - (n & (bpMu*bpGradDU.T() - (bpMu + bpLambda)*bpGradDU))
          - n*bpLambda*tr(bpGradDU)
        )
       /(2.0*bpMu + bpLambda);
    }
}
