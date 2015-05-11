{
    scalar CoNumSigma = max
    (
        sqrt
        (
            2*M_PI*sigma*sqr(aMesh.edgeInterpolation::deltaCoeffs())
            *aMesh.edgeInterpolation::deltaCoeffs()
            /rhol
        )
    ).value()*runTime.deltaT().value();

    Info<< "Max Capillary Courant Number = " << CoNumSigma << '\n' << endl;
}
