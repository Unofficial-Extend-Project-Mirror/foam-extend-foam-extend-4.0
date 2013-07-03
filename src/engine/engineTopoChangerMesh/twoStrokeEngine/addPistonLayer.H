if (piston().patchID().active())
{
    // Add faces for piston layering
    // Note: because of operation of layer addition/removal
    // (reduce function in layer addition/removal thickness)
    // the layering modifier needs to be present on all processors
    // even if the patch size is zero
    // HJ, 7/Mar/2011

    faceSet pistonFaceSet(*this, piston().pistonFaceSetName());

    boolList flipPistonFaces(pistonFaceSet.toc().size(), false);
    label nSet = pistonFaceSet.toc().size();
    label nFlip = 0;

    forAll (flipPistonFaces, facei)
    {
        scalar scalProd =
            (faceAreas()[pistonFaceSet.toc()[facei]] & vector(0, 0, 1));

        if (scalProd < 0)
        {
            flipPistonFaces[facei] = true;
            nFlip++;
        }
    }

    Info << "nSet = " << nSet << endl;
    Info << "nFlip = " << nFlip << endl;

    fz[nFaceZones] =
    (
        new faceZone
        (
            "pistonLayerFaces",
            pistonFaceSet.toc(),
            flipPistonFaces,
            nFaceZones,
            faceZones()
        )
    );

    nFaceZones++;

    cellSet movingPistonCells(*this, piston().pistonCellSetName());

    Info<< "Adding piston cell set" << endl;
    cz[nCellZones] =
    (
        new cellZone
        (
            "pistonCells",
            movingPistonCells.toc(),
            nCellZones,
            cellZones()
        )
    );

    nCellZones++;
}
