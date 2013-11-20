if (mesh.moving())
{
    dimensionedScalar one("one", dimless, 1.0);

    volScalarField volConservation
    (
        "volConservation",
        -fvc::div(mesh.phi())
    );

    // The ddt term constructed by hand because it must be Euler
    volConservation.internalField() +=
        (1.0 - mesh.V0()/mesh.V())/runTime.deltaT().value();

    scalar sumLocalVolContErr = runTime.deltaT().value()*
        mag(volConservation)().weightedAverage(mesh.V()).value();

    scalar globalVolContErr = runTime.deltaT().value()*
        volConservation.weightedAverage(mesh.V()).value();

    Info<< "volume continuity errors : "
        << "volume = " << sum(mesh.V()).value()
        << ", max error = " << max(volConservation.internalField())
        << ", sum local = " << sumLocalVolContErr
        << ", global = " << globalVolContErr << endl;
}
