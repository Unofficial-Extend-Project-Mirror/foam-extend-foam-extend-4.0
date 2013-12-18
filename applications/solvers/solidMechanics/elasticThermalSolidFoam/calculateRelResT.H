{
  // temperature relative residual
  scalar maxDT = gMax(mag(T.internalField() - T.oldTime().internalField()));
    relResT =
        gMax
        (
            mag(T.internalField() - T.prevIter().internalField())
           /(maxDT + SMALL)
        );
}
