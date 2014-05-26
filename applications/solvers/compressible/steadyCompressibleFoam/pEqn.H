{
    volScalarField rUA = 1.0/UEqn.A();

    for (int corr = 0; corr < nCorr; corr++)
    {
        U = rUA*UEqn.H();

        // Execute ddtPhiCorr before recalculating flux
        // HJ, 27/Apr/2010
        surfaceScalarField phid
        (
            "phid",
            fvc::interpolate(psi)*
            (
                (fvc::interpolate(U) & mesh.Sf())
              + fvc::ddtPhiCorr(rUA, rho, U, phi)
            )
        );

        // Calculate phi for boundary conditions
        phi = fvc::interpolate(rho*U) & mesh.Sf();

        p.storePrevIter();

        for (int nonOrth = 0; nonOrth <= nNonOrthCorr; nonOrth++)
        {
            fvScalarMatrix pEqn
            (
                fvm::ddt(psi, p)
              + fvm::div(phid, p)
              - fvm::laplacian(rho*rUA, p)
            );

            // Retain the residual from the first pressure solution
            eqnResidual = pEqn.solve().initialResidual();

            if (corr == 0 && nonOrth == 0)
            {
                maxResidual = max(eqnResidual, maxResidual);
            }

            // Calculate the flux
            if (nonOrth == nNonOrthCorr)
            {
                phi = pEqn.flux();
            }
        }

#       include "compressibleContinuityErrs.H"

        // Relax the pressure
        p.relax();

        U -= rUA*fvc::grad(p);
        U.correctBoundaryConditions();
    }

    // Bound the pressure
    if (min(p) < pMin || max(p) > pMax)
    {
        p.max(pMin);
        p.min(pMax);
        p.correctBoundaryConditions();
    }

    // Bound the velocity
    volScalarField magU = mag(U);

    if (max(magU) > UMax)
    {
        volScalarField Ulimiter = pos(magU - UMax)*UMax/(magU + smallU)
            + neg(magU - UMax);
        Ulimiter.max(scalar(0));
        Ulimiter.min(scalar(1));

        U *= Ulimiter;
        U.correctBoundaryConditions();
    }
}
