
{
    label patchID = mesh.boundaryMesh().findPatchID("freeSurface");

    if(patchID == -1)
    {
        FatalErrorIn(args.executable())
            << "Can't find patch with name freeSurface" << abort(FatalError);
    }

    if
    (
        DU.boundaryField()[patchID].type()
     != tractionDisplacementIncrementFvPatchVectorField::typeName
    )
    {
        FatalErrorIn(args.executable())
            << "Bounary condition on " << DU.name()
                <<  " is "
                << DU.boundaryField()[patchID].type()
                << "for patch: " << mesh.boundaryMesh()[patchID].name()
                <<  ", instead "
                << tractionDisplacementIncrementFvPatchVectorField::typeName
                << abort(FatalError);
    }

    tractionDisplacementIncrementFvPatchVectorField& tractionPatch =
        refCast<tractionDisplacementIncrementFvPatchVectorField>
        (
            DU.boundaryField()[patchID]
        );

    vectorField n = mesh.boundary()[patchID].nf();

    vector traction = vector::zero;
    vector DTraction = vector::zero;

    vectorField relaxedTraction =
        (
            n
          & (
                sigma.boundaryField()[patchID]
              + DSigmaCorr.boundaryField()[patchID]
            )
        );

    tractionPatch.DTraction() = DTraction + (traction - relaxedTraction);

    tractionPatch.DPressure() = 0.0;
}


