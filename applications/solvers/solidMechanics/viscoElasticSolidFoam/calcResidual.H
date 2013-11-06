{
    scalarField magDU = mag(DU.internalField());

    forAll(magDU, cellI)
    {
        if (magDU[cellI] < SMALL)
        {
            magDU[cellI] = SMALL;
        }
    }

    residual =
        gMax
        (
            mag
            (
                DU.internalField()
              - DU.prevIter().internalField()
            )
           /magDU
        );
}
