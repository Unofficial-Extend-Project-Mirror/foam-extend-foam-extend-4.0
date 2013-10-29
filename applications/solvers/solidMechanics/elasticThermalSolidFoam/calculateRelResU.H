{
    // displacement
    scalar maxDU = gMax(mag(U.internalField() - U.oldTime().internalField()));
    relResU =
        gMax
        (
            mag(U.internalField() - U.prevIter().internalField())
           /(maxDU + SMALL)
        );
}
