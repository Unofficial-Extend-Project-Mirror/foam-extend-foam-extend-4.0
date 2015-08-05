    scalar DiNum = -GREAT;

    forAll(solidRegions, regionI)
    {
        DiNum = max
        (
            solidRegionDiffNo
            (
                solidRegions[regionI],
                runTime,
                rhosCps[regionI],
                Kappas[regionI]
            ),
            DiNum
        );
    }
