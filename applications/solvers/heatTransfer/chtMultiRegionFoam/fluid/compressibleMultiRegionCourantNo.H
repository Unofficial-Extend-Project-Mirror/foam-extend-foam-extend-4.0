    scalar CoNum = -GREAT;
    forAll(fluidRegions, regionI)
    {
        CoNum = max
        (
            compressibleCourantNo
            (
                fluidRegions[regionI],
                runTime,
                rhof[regionI],
                phif[regionI]
            ),
            CoNum
        );
    }
