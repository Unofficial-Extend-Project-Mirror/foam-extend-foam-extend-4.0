{
  scalar magDU = gMax(mag(DU.internalField()));

  if (magDU < SMALL)
    {
      magDU = SMALL;
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
