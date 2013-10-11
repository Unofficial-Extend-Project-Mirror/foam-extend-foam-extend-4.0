{
    scalar maxDU = gMax(mag(DU.internalField()));

    relativeResidual =
        gMax
        (
            mag(DU.internalField() - DU.prevIter().internalField())
           /(maxDU + SMALL)
        );

    if (lduMatrix::debug == 1)
    {
        Info << "Relative residual = " << relativeResidual << endl;
    }
}
