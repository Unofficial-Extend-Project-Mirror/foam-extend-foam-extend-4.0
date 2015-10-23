if (dynamicTimeStep)
{
    if
    (
        //(maxEffTraction < 0.999*CzLaw.sigmaMax().value())
        (returnReduce(maxEffTractionFraction, maxOp<scalar>()) < 0.99)
      //&& (cohesivePatchU.size() == 0)
      && (mag(runTime.deltaT().value() - deltaTmax) < SMALL)
    )
    {
        runTime.setDeltaT(deltaTmax);
    }
    else
    {
        scalar newDeltaT = deltaTmin;

        if (newDeltaT/runTime.deltaT().value() < 0.5)
        {
            newDeltaT = 0.5*runTime.deltaT().value();
            Info << "Reducing time step" << nl;
        }

        runTime.setDeltaT(newDeltaT);
    }

    Pout << "Current time step size: "
        << runTime.deltaT().value() << " s" << endl;

    scalar maxDT = runTime.deltaT().value();

    if(mag(returnReduce(maxDT, maxOp<scalar>()) - runTime.deltaT().value()) > SMALL)
    {
        FatalError << "Processors have different time-steps!"
            << exit(FatalError);
    }
}
