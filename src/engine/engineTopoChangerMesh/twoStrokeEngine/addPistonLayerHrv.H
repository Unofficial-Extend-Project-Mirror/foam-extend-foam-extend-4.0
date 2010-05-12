if (piston().patchID().active())
{
    // add faces for piston layering

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

    fz.append
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

//     {
//         pointSet movingPistonPoints(*this, piston().pistonPointSetName());

//         pz.append
//         (
//             new pointZone
//             (
//                 "pistonPoints",
//                 movingPistonPoints.toc(),
//                 nPointZones,
//                 pointZones()
//             )
//         );

//         nPointZones++;
//     }

    cellSet movingPistonCells(*this, piston().pistonCellSetName());

    Info<< "Adding piston cell set" << endl;
    cz.append
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
