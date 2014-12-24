
{
    label patchID = mesh.boundaryMesh().findPatchID("loadedEnd");

    if(patchID == -1)
    {
        FatalErrorIn(args.executable())
            << "Can't find patch with name loadedEnd" << abort(FatalError);
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

    if(runTime.timeIndex()==1)
    {
        traction = vector::zero;
        DTraction = maxTraction;
    }
    else if(runTime.timeIndex()>=2)
    {
        traction = maxTraction;
        DTraction = vector::zero;
    }

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


