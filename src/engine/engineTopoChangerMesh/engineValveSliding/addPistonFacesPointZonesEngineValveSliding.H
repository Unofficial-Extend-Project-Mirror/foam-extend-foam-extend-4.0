// Add the piston points zone
if (piston().patchID().active())
{
/*
    pointSet movingPistonPoints(*this, piston().pistonPointSetName());

    pz.append
    (
        new pointZone
        (
            "movingPistonPoints",
            movingPistonPoints.toc(),
            nPointZones,
            pointZones()
        )
    );
    nPointZones++;
*/

    cellSet movingPistonCells(*this, piston().pistonCellSetName());

    cz.append
    (
        new cellZone
        (
            "movingPistonCells",
            movingPistonCells.toc(),
            nCellZones,
            cellZones()
        )
    );

    nCellZones++;
}
