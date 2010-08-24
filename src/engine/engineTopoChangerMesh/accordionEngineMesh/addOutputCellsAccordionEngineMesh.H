{
    word cylinderOutputName
    (
        engTime().engineDict().lookup("cylinderSetName")
    );

    cellSet outputCellSet(*this, cylinderOutputName);

    cz.append
    (
        new cellZone
        (
            "cylinderOutputCells",
            outputCellSet.toc(),
            nCellZones,
            cellZones()
        )
    );

    nCellZones++;
}
