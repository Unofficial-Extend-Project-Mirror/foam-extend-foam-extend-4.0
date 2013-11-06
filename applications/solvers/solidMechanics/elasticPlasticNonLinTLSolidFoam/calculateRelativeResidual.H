{
  scalar maxDU = gMax(mag(DU.internalField()));

    relativeResidual =
        gMax
      (
       mag(DU.internalField() - DU.prevIter().internalField())
       /(maxDU + SMALL)
       );

    //resid = (DU - DU.prevIter())/(mag(DU)+dimensionedScalar("small",dimLength,SMALL));
}
