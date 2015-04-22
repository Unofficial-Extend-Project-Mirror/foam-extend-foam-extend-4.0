{
    // Bound the velocity
    volScalarField magU = mag(U);

    if (max(magU) > UMax)
    {
        volScalarField Ulimiter =
            pos(magU - UMax)*UMax/(magU + smallU)
          + neg(magU - UMax);

        Ulimiter.max(scalar(0));
        Ulimiter.min(scalar(1));

        U *= Ulimiter;
        U.correctBoundaryConditions();
    }
}
