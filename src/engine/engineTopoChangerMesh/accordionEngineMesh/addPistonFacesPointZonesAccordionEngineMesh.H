
    // Add the piston points zone
    if (piston().patchID().active())
    {
        pointSet movingPistonPoints(*this, pistonAuxPoints_);

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
    }
