    scalar newTotalVolume = gSum(mesh.cellVolumes());
    scalar totalVolRatio = newTotalVolume/totalVolume;

    Info<< "Volume: new = " << newTotalVolume << " old = " << totalVolume
        << " change = " << Foam::mag(newTotalVolume - totalVolume)
        << " ratio = " << totalVolRatio - 1 << endl;

    totalVolume = newTotalVolume;
