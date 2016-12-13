{
    U = UEqn.H()/UEqn.A();

#   include "limitU.H"

    while (pimple.correctNonOrthogonal())
    {
        // Calculate phi for boundary conditions
        phi = rhof*
        (
            (fvc::interpolate(U) & mesh.Sf())
          - fvc::meshPhi(rho, U)
        );

        surfaceScalarField phid2 = rhoReff/rhof*phi;

        surfaceScalarField phid("phid", psisf/rhof*phi);

        // Store pressure for under-relaxation
        p.storePrevIter();

        volScalarField divPhid
        (
            "divPhid",
            fvc::div(phid)
        );

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

        pEqn.solve
        (
            mesh.solutionDict().solver(p.select(pimple.finalInnerIter()))
        );

        // Calculate the flux
        if (pimple.finalNonOrthogonalIter())
        {
            phi = phid2 + pEqn.flux();
        }


        // Bound the pressure
        if (min(p) < pMin || max(p) > pMax)
        {
            p.max(pMin);
            p.min(pMax);
            p.correctBoundaryConditions();
        }

        // Relax the pressure
        p.relax();
    }

#   include "compressibleContinuityErrs.H"

    U -= fvc::grad(p)/UEqn.A();
    U.correctBoundaryConditions();

#   include "limitU.H"
}
