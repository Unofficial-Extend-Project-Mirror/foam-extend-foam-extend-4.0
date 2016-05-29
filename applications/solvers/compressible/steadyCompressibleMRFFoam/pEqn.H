{
    volScalarField rUA = 1.0/UEqn.A();

    surfaceScalarField psisf = fvc::interpolate(psis);
    surfaceScalarField rhof = fvc::interpolate(rho);

    // Needs to be outside of loop since p is changing, but psi and rho are not.
    surfaceScalarField rhoReff = rhof - psisf*fvc::interpolate(p);

    while (pimple.correct())
    {
        U = rUA*UEqn.H();

        // Calculate phi for boundary conditions
        phi = rhof*fvc::interpolate(U) & mesh.Sf();

        surfaceScalarField phid2 = rhoReff/rhof*phi;

        surfaceScalarField phid("phid", psisf/rhof*phi);

        // Make fluxes relative within the MRF zone
        mrfZones.relativeFlux(rhoReff, phi);
        mrfZones.relativeFlux(psisf, phid);
        mrfZones.relativeFlux(rhoReff, phid2);

        p.storePrevIter();

        while (pimple.correctNonOrthogonal())
        {
            fvScalarMatrix pEqn
            (
                fvm::ddt(psis, p)
              + fvm::div(phid, p)
                // Convective flux relaxation terms
              + fvm::SuSp(-divPhid, p)
              + divPhid*p
              + fvc::div(phid2)
              - fvm::laplacian(rho*rUA, p)
            );

            pEqn.solve();

            // Calculate the flux
            if (pimple.finalNonOrthogonalIter())
            {
                phi = phid2 + pEqn.flux();
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
