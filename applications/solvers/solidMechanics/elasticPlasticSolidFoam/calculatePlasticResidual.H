{
  scalar magDEp = gMax(mag(DEpsilonP.internalField()));

  if (magDEp < SMALL)
    {
      magDEp = SMALL;
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
