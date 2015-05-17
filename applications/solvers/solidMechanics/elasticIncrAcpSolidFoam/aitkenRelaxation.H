// aitken acceleration
aitkenDelta.storePrevIter();

// update delta
aitkenDelta = (DU - DU.prevIter()) / aitkenInitialRes;

// update relaxation factor
if(iCorr == 0)
{
    aitkenTheta = 0.01;
}
else
{
    vectorField b = aitkenDelta.internalField() - aitkenDelta.prevIter().internalField();
    //scalar sumMagB = gSum(mag(b));
    scalar sumMagB = gSum(magSqr(b));
    if(sumMagB < SMALL)
    {
        //Warning << "Aitken under-relaxation: denominator less then SMALL"
        //    << endl;
        sumMagB += SMALL;
    }

    aitkenTheta = -aitkenTheta*
        gSum(aitkenDelta.prevIter().internalField() & b)/sumMagB;
}

// correction to the latest DU
DU += aitkenTheta*aitkenDelta*aitkenInitialRes;
