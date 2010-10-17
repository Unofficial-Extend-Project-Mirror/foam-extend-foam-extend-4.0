{
    // Detach patches
#   include "detachPatches.H"

#   include "readPISOControls.H"
#   include "CourantNo.H"

    fvVectorMatrix UEqn
    (
        fvm::ddt(U)
      + fvm::div(phi, U)
      - fvm::laplacian(nu, U)
    );

    solve(UEqn == -fvc::grad(p));

    // --- PISO loop

    for (int corr = 0; corr < nCorr; corr++)
    {
        U = UEqn.H()/UEqn.A();
        U.correctBoundaryConditions();

        phi = fvc::interpolate(U) & mesh.Sf();

        for (int nonOrth = 0; nonOrth <= nNonOrthCorr; nonOrth++)
        {
            fvScalarMatrix pEqn
            (
                fvm::laplacian(1.0/UEqn.A(), p) == fvc::div(phi)
            );

            pEqn.setReference(pRefCell, pRefValue);
            pEqn.solve();

            if (nonOrth == nNonOrthCorr)
            {
                phi -= pEqn.flux();
            }
        }

#       include "continuityErrs.H"

        U -= fvc::grad(p)/UEqn.A();
        U.correctBoundaryConditions();
    }
}
