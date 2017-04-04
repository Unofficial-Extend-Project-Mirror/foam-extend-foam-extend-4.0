{
    volScalarField rUrelA = 1.0/UrelEqn.A();

    surfaceScalarField psisf = fvc::interpolate(psis);
    surfaceScalarField rhof = fvc::interpolate(rho);

    // Needs to be outside of loop since p is changing, but psi and rho are not.
    surfaceScalarField rhoReff = rhof - psisf*fvc::interpolate(p);

    while (pimple.correct())
    {
        Urel = rUrelA*UrelEqn.H();

        // Calculate phi for boundary conditions
        phi = rhof*fvc::interpolate(Urel) & mesh.Sf();

        surfaceScalarField phid2 = rhoReff/rhof*phi;

        surfaceScalarField phid("phid", psisf/rhof*phi);

        p.storePrevIter();

        while (pimple.correctNonOrthogonal())
        {
            fvScalarMatrix pEqn
            (
                fvm::div(phid, p)
              + fvc::div(phid2)
              - fvm::laplacian(rho*rUrelA, p)
            );

            pEqn.solve();

            // Calculate the flux
            if (pimple.finalNonOrthogonalIter())
            {
                phi = phid2 + pEqn.flux();
            }
        }

        // Use incompressible continuity error check: div(rho U) = 0
#       include "continuityErrs.H"

        // Relax the pressure
        p.relax();

        Urel -= rUrelA*fvc::grad(p);
        Urel.correctBoundaryConditions();
    }

    // Bound the pressure
    if (min(p) < pMin || max(p) > pMax)
    {
        p.max(pMin);
        p.min(pMax);
        p.correctBoundaryConditions();
    }

    // Bound the velocity
    volScalarField magUrel = mag(Urel);

    if (max(magUrel) > UrelMax)
    {
        volScalarField Urellimiter =
            pos(magUrel - UrelMax)*UrelMax/(magUrel + smallUrel)
            + neg(magUrel - UrelMax);
        Urellimiter.max(scalar(0));
        Urellimiter.min(scalar(1));

        Urel *= Urellimiter;
        Urel.correctBoundaryConditions();
    }
}
