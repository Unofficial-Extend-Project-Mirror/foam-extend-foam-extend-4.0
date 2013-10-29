{
    scalarField magDU = mag(DU.internalField());

    forAll(magDU, cellI)
    {
        if (magDU[cellI] < SMALL)
        {
            magDU[cellI] = SMALL;
        }
    }

    relativeResidual =
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
