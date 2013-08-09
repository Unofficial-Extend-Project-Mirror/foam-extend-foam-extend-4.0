    fvMesh solidMesh
    (
        IOobject
        (
            "solid",
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ
        )
    );

{
    const polyPatchList& patches = solidMesh.boundaryMesh();

    forAll (patches, patchI)
    {
        if (isType<regionCouplePolyPatch>(patches[patchI]))
        {
            const regionCouplePolyPatch& rcp =
                refCast<const regionCouplePolyPatch>(patches[patchI]);

            // Attach it here
            rcp.attach();
        }
    }

    // Force recalculation of weights
    solidMesh.surfaceInterpolation::movePoints(); 
}

