{
    scalar t = runTime.value();
    scalar tNext = t + runTime.deltaT().value();

    instantList Times = runTime.times();

    sigma = dimensionedSymmTensor("zero", dimForce/dimArea, symmTensor::zero);
    DSigmaCorr =
        dimensionedSymmTensor
        (
            "zero",
            dimForce/dimArea,
            symmTensor::zero
        );

    for (label i=1; i<Times.size(); i++)
    {
        runTime.setTime(Times[i], i);

        if(runTime.timeIndex() != i)
        {
            FatalErrorIn(args.executable())
                << "The strain increment DEpsilon must be stored for "
                << "each calculated time step. "
                << exit(FatalError);
        }

        IOobject DEpsilonHeader
        (
            "DEpsilon",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );

        // Check DEpsilon exists
        if (DEpsilonHeader.headerOk())
        {
            volSymmTensorField DEpsilonOld(DEpsilonHeader, mesh);

            scalar tau = runTime.value() - m*runTime.deltaT().value();

            if(tau < 0)
            {
                sigma += 2.0*rheology.mu(t)*DEpsilonOld
                    + rheology.lambda(t)*(I*tr(DEpsilonOld));

                DSigmaCorr += 2.0*rheology.mu(tNext)*DEpsilonOld
                    + rheology.lambda(tNext)*(I*tr(DEpsilonOld));
            }
            else
            {
                sigma += 2.0*rheology.mu(t - tau)*DEpsilonOld
                    + rheology.lambda(t - tau)*(I*tr(DEpsilonOld));

                DSigmaCorr += 2.0*rheology.mu(tNext - tau)*DEpsilonOld
                    + rheology.lambda(tNext - tau)*(I*tr(DEpsilonOld));
            }
        }
        else
        {
            Info << "No DEpsilon" << endl;
        }
    }

    if(Times.size()>=2)
    {
        runTime++;
    }

    scalar tau = runTime.value() - m*runTime.deltaT().value();

    sigma += 2.0*rheology.mu(t - tau)*DEpsilon
        + rheology.lambda(t - tau)*(I*tr(DEpsilon));

    DSigmaCorr += 2.0*rheology.mu(tNext - tau)*DEpsilon
        + rheology.lambda(tNext - tau)*(I*tr(DEpsilon));

    DSigmaCorr -= sigma;
}
