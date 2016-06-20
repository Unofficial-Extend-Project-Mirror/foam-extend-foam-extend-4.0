if(fsi.closedFluidDomain())
{
    scalar inletFlux = 0;
    forAll(U_.boundaryField(), patchI)
    {
        if (patchI != fsi.fluidPatchIndex())
        {
            inletFlux +=
                gSum
                (
                    U_.boundaryField()[patchI]
                  & mesh.Sf().boundaryField()[patchI]
                );
        }
    }

    phi_.boundaryField()[fsi.fluidPatchIndex()] =
    (
        U_.boundaryField()[fsi.fluidPatchIndex()]
      & mesh.Sf().boundaryField()[fsi.fluidPatchIndex()]
    );

    Info << "Inlet flux: " << inletFlux << endl;

    Info << "Global continuity error (before): "
        << gSum(phi_.boundaryField()[fsi.fluidPatchIndex()]) + inletFlux
        << endl;

    scalarField weights =
        mag(phi_.boundaryField()[fsi.fluidPatchIndex()]);

    if(mag(gSum(weights)) > VSMALL)
    {
        weights /= gSum(weights) + SMALL;
    }

    phi_.boundaryField()[fsi.fluidPatchIndex()] -=
        weights
       *(gSum(phi_.boundaryField()[fsi.fluidPatchIndex()]) + inletFlux);

    Info << "Global continuity error (after): "
        << gSum(phi_.boundaryField()[fsi.fluidPatchIndex()]) + inletFlux
        << endl;

    // Correct velocity

    const vectorField n = mesh.boundary()[fsi.fluidPatchIndex()].nf();

    U_.boundaryField()[fsi.fluidPatchIndex()] ==
        U_.boundaryField()[fsi.fluidPatchIndex()]
      - n*(n & U_.boundaryField()[fsi.fluidPatchIndex()])
      + n
       *(
           phi_.boundaryField()[fsi.fluidPatchIndex()]
          /mesh.magSf().boundaryField()[fsi.fluidPatchIndex()]
        );
}
