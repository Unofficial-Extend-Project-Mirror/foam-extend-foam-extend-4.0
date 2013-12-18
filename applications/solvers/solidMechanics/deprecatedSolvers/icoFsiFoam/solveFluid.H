{
    // SIMPLE loop

    for (int corr=0; corr<nCorr; corr++)
    {
        fvVectorMatrix UEqn
        (
            fvm::ddt(U)
          + fvm::div(phi, U)
          - fvm::laplacian(nu, U)
        );

        UEqn.relax();

        solve(UEqn == -fvc::grad(p));

        U = UEqn.H()/UEqn.A();
        U.correctBoundaryConditions();

        adjustPhi(phi, U, p);

        phi = fvc::interpolate(U) & mesh.Sf();

        p.storePrevIter();

        for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
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

        p.relax();

#       include "movingMeshContinuityErrs.H"

        U -= fvc::grad(p)/UEqn.A();
        U.correctBoundaryConditions();

        // Make the fluxes relative
        phi -= fvc::meshPhi(U);
    }
}
