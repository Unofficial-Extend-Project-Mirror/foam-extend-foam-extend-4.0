{
    // Bound the pressure
    dimensionedScalar p1 = min(p);
    dimensionedScalar p2 = max(p);

    if (p1 < pMin || p2 > pMax)
    {
        Info<< "p: " << p1.value() << " " << p2.value()
            << ".  Bounding." << endl;

        p.max(pMin);
        p.min(pMax);
        p.correctBoundaryConditions();
    }

    // Bound the velocity
    volScalarField magU = mag(U);
    dimensionedScalar U1 = max(magU);

    if (U1 > UMax)
    {
        Info<< "U: " << U1.value() << ".  Bounding." << endl;

        volScalarField Ulimiter = pos(magU - UMax)*UMax/(magU + smallU)
            + neg(magU - UMax);
        Ulimiter.max(scalar(0));
        Ulimiter.min(scalar(1));

        U *= Ulimiter;
        U.correctBoundaryConditions();
    }
}
