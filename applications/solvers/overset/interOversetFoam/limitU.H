{
    scalar limitMagU = readScalar(pimple.dict().lookup("limitMagU"));

    volScalarField magU(mag(U));

    scalar maxMagU = max(magU).value();

    Info<< "mag(U): max: " << maxMagU
        << " avg: " << magU.weightedAverage(mesh.V()).value();

    if (maxMagU > limitMagU)
    {
        U.internalField() *=
            neg(magU.internalField() - limitMagU)
          + pos(magU.internalField() - limitMagU)*
            limitMagU/(magU.internalField() + SMALL);
        U.correctBoundaryConditions();
        oversetFvPatchVectorField::oversetInterpolate(U); // Overset update
        // Note: if implicit fringe is switched on, we are doing the
        // interpolation twice (once in correctBoundaryConditions and once in
        // oversetInterpolate). Reorganize. VV, 4/Oct/2016.

        Info << " ...limiting" << endl;
    }
    else
    {
        Info << endl;
    }
}
