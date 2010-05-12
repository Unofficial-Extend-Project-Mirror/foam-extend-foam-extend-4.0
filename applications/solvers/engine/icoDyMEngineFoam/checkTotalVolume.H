{
    // Volume and motion flux check

    scalar newTotalVolume = sum(mesh.V()).value();

    volScalarField motionContErr =
        fvc::ddt(dimensionedScalar("1", dimless, 1.0), mesh)
      - fvc::div(fvc::meshPhi(U));

    scalar sumLocalContErr = runTime.deltaT().value()*
        mag(motionContErr)().weightedAverage(mesh.V()).value();

    scalar maxContErr =
        runTime.deltaT().value()*max(motionContErr.internalField());

    Info<< "Volume: new = " << newTotalVolume << " old = " << totalVolume
        << " change = " << Foam::mag(newTotalVolume - totalVolume) << nl
        << "Motion continuity errors : sum local = " << sumLocalContErr
        << ", maximum = " << maxContErr
        << endl;

    totalVolume = newTotalVolume;
}
