label pRefCell = 0;
scalar pRefValue = 0.0;
bool pNeedRef = false;
bool procHasRef = false;

if (interface.twoFluids())
{
    pNeedRef = true;
}

forAll (p.boundaryField(), patchI)
{
    // Check only fluid B
    if (fluidIndicator.boundaryField()[patchI].size())
    {
        if (min(fluidIndicator.boundaryField()[patchI]) < SMALL)
        {
            if (p.boundaryField()[patchI].fixesValue())
            {
                pNeedRef = false;
            }
        }
    }
}

pNeedRef = returnReduce<bool>(pNeedRef, andOp<bool>());

// Find reference cell
if (pNeedRef)
{
    point refPointi(mesh.solutionDict().subDict("PISO").lookup("pRefPoint"));
    label refCelli = mesh.findCell(refPointi);
    label hasRef = (refCelli >= 0 ? 1 : 0);
    label sumHasRef = returnReduce<label>(hasRef, sumOp<label>());

    if (sumHasRef != 1)
    {
        FatalErrorIn(args.executable())
            << "Unable to set reference cell for field "
                << p.name()
                << nl << "    Reference point pRefPoint"
                << " found on " << sumHasRef << " domains (should be one)"
                << nl << exit(FatalError);
    }

    if (hasRef)
    {
        pRefCell = refCelli;
        procHasRef = true;
    }

    pRefValue = 
        readScalar(mesh.solutionDict().subDict("PISO").lookup("pRefValue"));
}

