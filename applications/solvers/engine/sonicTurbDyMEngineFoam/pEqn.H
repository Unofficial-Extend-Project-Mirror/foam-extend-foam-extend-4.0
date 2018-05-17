{
    rUA = 1.0/UEqn.A();
    U = rUA*UEqn.H();

    if (pimple.transonic())
    {
        surfaceScalarField phid =
            fvc::interpolate(thermo.psi())*
            ((fvc::interpolate(U) & mesh.Sf()) - fvc::meshPhi(rho, U));

        while (pimple.correctNonOrthogonal())
        {
            fvScalarMatrix pEqn
            (
                fvm::ddt(psi, p)
              + fvm::div(phid, p, "div(phid,p)")
              - fvm::laplacian(rho*rUA, p)
            );

            pEqn.solve();

            if (pimple.finalNonOrthogonalIter())
            {
                phi == pEqn.flux();
            }
        }
    }
    else
    {
        phi = fvc::interpolate(rho)*
            ((fvc::interpolate(U) & mesh.Sf()) - fvc::meshPhi(rho, U));

        while (pimple.correctNonOrthogonal())
        {
            fvScalarMatrix pEqn
            (
                fvm::ddt(psi, p)
              + fvc::div(phi)
              - fvm::laplacian(rho*rUA, p)
            );

            pEqn.solve();

            if (pimple.finalNonOrthogonalIter())
            {
                phi += pEqn.flux();
            }
        }
    }

    // Explicitly relax pressure except for last corrector
    if (!pimple.finalIter())
    {
        p.relax();
    }

#   include "compressibleContinuityErrs.H"

    U -= rUA*fvc::grad(p);
    U.correctBoundaryConditions();

    DpDt = fvc::DDt(surfaceScalarField("phiU", phi/fvc::interpolate(rho)), p);
}
