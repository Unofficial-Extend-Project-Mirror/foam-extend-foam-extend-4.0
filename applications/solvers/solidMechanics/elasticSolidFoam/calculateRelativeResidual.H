{
  //scalarField magDU = mag(U.internalField() - U.oldTime().internalField());
    // forAll(magDU, cellI)
    // {
    //     if (magDU[cellI] < SMALL)
    //     {
    //         magDU[cellI] = SMALL;
    //     }
    // }

    // relativeResidual =
    //     gMax
    //     (
    //         mag
    //         (
    //             U.internalField()
    //           - U.prevIter().internalField()
    //         )
    //        /magDU
    //     );

    scalar maxDU = gMax(mag(U.internalField() - U.oldTime().internalField()));
    relativeResidual =
        gMax
        (
            mag(U.internalField() - U.prevIter().internalField())
           /(maxDU + SMALL)
        );
}
