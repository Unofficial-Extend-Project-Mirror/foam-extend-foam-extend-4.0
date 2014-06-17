{
    Info << "\nCalculate total stress" << endl;

    scalar t = runTime.value();
    scalar tNext = t + runTime.deltaT().value();

    instantList Times = runTime.times();

//     Info << "Number of times: " << Times.size() << endl;

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
                << "The strain increment Depsilon must be stored for "
                    << "each calculated time step. "
                    << exit(FatalError);
        }

        IOobject DepsilonHeader
        (
            "Depsilon",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );

        // Check Depsilon exists
        if (DepsilonHeader.headerOk())
        {
            volSymmTensorField DepsilonOld(DepsilonHeader, mesh);

            scalar tau = runTime.value() - m*runTime.deltaT().value();

            if(tau < 0)
            {
                sigma += 2.0*rheology.mu(t)*DepsilonOld
                    + rheology.lambda(t)*(I*tr(DepsilonOld));

                DSigmaCorr += 2.0*rheology.mu(tNext)*DepsilonOld
                    + rheology.lambda(tNext)*(I*tr(DepsilonOld));
            }
            else
            {
                sigma += 2.0*rheology.mu(t - tau)*DepsilonOld
                    + rheology.lambda(t - tau)*(I*tr(DepsilonOld));

                DSigmaCorr += 2.0*rheology.mu(tNext - tau)*DepsilonOld
                    + rheology.lambda(tNext - tau)*(I*tr(DepsilonOld));
            }
        }
        else
        {
            Info << "No Depsilon" << endl;
        }
    }

    if(Times.size()>=2)
    {
        runTime++;
    }

    scalar tau = runTime.value() - m*runTime.deltaT().value();

    sigma += 2.0*rheology.mu(t - tau)*Depsilon
        + rheology.lambda(t - tau)*(I*tr(Depsilon));

    DSigmaCorr += 2.0*rheology.mu(tNext - tau)*Depsilon
        + rheology.lambda(tNext - tau)*(I*tr(Depsilon));

    DSigmaCorr -= sigma;

    //Info << "Current time = " << runTime.timeName() << endl;
}

