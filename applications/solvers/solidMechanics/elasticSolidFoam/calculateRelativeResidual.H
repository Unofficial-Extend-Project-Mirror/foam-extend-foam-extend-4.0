{
    scalarField magU = mag(U.internalField() - U.oldTime().internalField());

    forAll(magU, cellI)
    {
        if (magU[cellI] < SMALL)
        {
            magU[cellI] = SMALL;
        }
    }

    relativeResidual =
        gMax
        (
            mag
            (
                U.internalField()
              - U.prevIter().internalField()
            )
           /magU
        );
}
