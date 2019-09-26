{
    p.boundaryField().updateCoeffs();

    volScalarField rAU("rAU", 1.0/UEqn().A());
    surfaceScalarField rAUf("(1|A(U))", fvc::interpolate(rAU));

    U = rAU*UEqn().H();
    UEqn.clear();

    simple.calcSteadyConsistentFlux(phi, U);
    adjustPhi(phi, U, p);

    surfaceScalarField buoyancyPhi =
        rAUf*fvc::interpolate(rhok)*(g & mesh.Sf());
    phi += buoyancyPhi;

    while (simple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian
            (
                rAUf/simple.aCoeff(U.name()),
                p,
                "laplacian(rAU," + p.name() + ')'
            )
         ==
            fvc::div(phi)
        );

        pEqn.setReference(pRefCell, pRefValue);

        pEqn.solve();

        if (simple.finalNonOrthogonalIter())
        {
            // Calculate the conservative fluxes
            phi -= pEqn.flux();
        }

        // U += rAU*fvc::reconstruct((buoyancyPhi - pEqn.flux())/rAUf);
        // U.correctBoundaryConditions();
    }

#   include "continuityErrs.H"

    // Explicitly relax pressure for momentum corrector
    p.relax();

    // Correct the momentum source with the pressure gradient flux
    // calculated from the relaxed pressure
    U += rAU*(fvc::reconstruct(buoyancyPhi/rAUf) - fvc::grad(p));
    U.correctBoundaryConditions();
}
