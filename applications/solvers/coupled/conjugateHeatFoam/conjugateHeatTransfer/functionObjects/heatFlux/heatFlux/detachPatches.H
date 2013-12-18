{
    const polyPatchList& patches = mesh.boundaryMesh();

    forAll (patches, patchI)
    {
        if (isType<regionCouplePolyPatch>(patches[patchI]))
        {
            const regionCouplePolyPatch& rcp =
                refCast<const regionCouplePolyPatch>(patches[patchI]);

            // Detach it here
            rcp.detach();
        }
    }

    // Force recalculation of weights
    mesh.surfaceInterpolation::movePoints(); 
}
