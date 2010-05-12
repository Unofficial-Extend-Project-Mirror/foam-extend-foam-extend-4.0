    scalar newTotalVolume = sum(mesh.cellVolumes());

    Info<< "Volume: new = " << newTotalVolume << " old = " << totalVolume
        << " change = " << Foam::mag(newTotalVolume - totalVolume) << endl;

    totalVolume = newTotalVolume;
