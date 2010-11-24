// if (pNeedRef && Pstream::master())
if (pNeedRef && procHasRef)
{
    pEqn.source()[pRefCell] += 
        pEqn.diag()[pRefCell]*pRefValue;

    pEqn.diag()[pRefCell] +=
        pEqn.diag()[pRefCell];
}
