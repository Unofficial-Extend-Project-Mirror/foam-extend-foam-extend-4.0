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

// plastic residual
{
  scalarField magDEp = mag(DEpsilonP.internalField());

    forAll(magDEp, cellI)
    {
        if (magDEp[cellI] < SMALL)
        {
            magDEp[cellI] = SMALL;
        }
    }

    plasticResidual =
        gMax
        (
            mag
            (
                DEpsilonP.internalField()
              - DEpsilonP.prevIter().internalField()
            )
           /magDEp
        );
}
