    if (adjustTimeStep)
    {
        if (CoNum > SMALL)
        {
            runTime.setDeltaT
            (
                min
                (
                    maxCo*runTime.deltaT().value()/CoNum,
                    maxDeltaT
                )
            );
        }
    }

