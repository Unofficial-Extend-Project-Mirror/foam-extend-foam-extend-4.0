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
     == fixedValueFvPatchVectorField::typeName
    )
    {
        if (runTime.timeIndex() == 1)
        {
            DU.boundaryField()[patchID] == vector(0, -0.00025, 0);
        }
        else
        {
            DU.boundaryField()[patchID] == vector(0, 0, 0);
        }
    }
    else if
    (
        DU.boundaryField()[patchID].type()
     == directionMixedFvPatchVectorField::typeName
    )
    {
        directionMixedFvPatchVectorField& topPatchU =
            refCast<directionMixedFvPatchVectorField>
            (
                DU.boundaryField()[patchID]
            );

        vectorField n = mesh.boundary()[patchID].nf();

        topPatchU.valueFraction() = sqr(n);

        if (runTime.timeIndex() == 1)
        {
            topPatchU.refValue() = vector(0, -0.00025, 0);
        }
        else
        {
            topPatchU.refValue() = vector(0, 0, 0);
        }
    }
    else
    {
        FatalErrorIn(args.executable())
            << "Bounary condition on " << DU.name()
                <<  " is "
                << DU.boundaryField()[patchID].type()
                << "for patch: " << mesh.boundaryMesh()[patchID].name()
                <<  ", instead "
                << fixedValueFvPatchVectorField::typeName
                << abort(FatalError);
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

        bottomPatchU.valueFraction() = sqr(n);

        bottomPatchU.refValue() = vector(0, 0, 0);
    }
}
