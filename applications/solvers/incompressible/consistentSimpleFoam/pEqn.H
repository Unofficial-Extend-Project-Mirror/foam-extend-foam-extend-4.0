{
    // Prepare clean 1/Ap without contribution from under-relaxation
    // HJ, 26/Oct/2015
    volScalarField rUA
    (
        "(1|A(U))",
        1/HUEqn().A()
    );

    // Store velocity under-relaxation point before using U for
    // the flux precursor
    U.storePrevIter();

    U = rUA*HUEqn().H();
    HUEqn.clear();
    phi = fvc::interpolate(U) & mesh.Sf();

    // Global flux balance
    adjustPhi(phi, U, p);

    // Non-orthogonal pressure corrector loop
    for (int nonOrth = 0; nonOrth <= nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(rUA, p) == fvc::div(phi)
        );

        pEqn.setReference(pRefCell, pRefValue);

        // Retain the residual from the first iteration
        if (nonOrth == 0)
        {
            eqnResidual = pEqn.solve().initialResidual();
            maxResidual = max(eqnResidual, maxResidual);
        }
        else
        {
            pEqn.solve();
        }

        if (nonOrth == nNonOrthCorr)
        {
            phi -= pEqn.flux();
        }
    }

#   include "continuityErrs.H"

    // Explicitly relax pressure for momentum corrector
    p.relax();

    // Momentum corrector
    // Note: since under-relaxation does not change aU, H/a in U can be
    // re-used.  HJ, 22/Jan/2016
    U = UUrf*(U - rUA*fvc::grad(p)) + (1 - UUrf)*U.prevIter();
    U.correctBoundaryConditions();
}
