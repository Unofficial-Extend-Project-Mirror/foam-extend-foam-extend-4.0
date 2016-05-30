    p.boundaryField().updateCoeffs();

    // Prepare clean 1/Ap without contribution from under-relaxation
    // HJ, 26/Oct/2015
    volScalarField rUA
    (
        "(1|A(Urel))",
        1/HUrelEqn().A()
    );

    // Store velocity under-relaxation point before using U for the flux
    // precursor
    Urel.storePrevIter();

    Urel = rUA*HUrelEqn().H();
    HUrelEqn.clear();
    phi = fvc::interpolate(Urel) & mesh.Sf();

    // Global flux balance
    adjustPhi(phi, Urel, p);

    // Non-orthogonal pressure corrector loop
    for (int nonOrth = 0; nonOrth <= nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(rUA, p) == fvc::div(phi)
        );

        pEqn.setReference(pRefCell, pRefValue);
        pEqn.solve();

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
    Urel = UUrf*(Urel - rUA*fvc::grad(p)) + (1 - UUrf)*Urel.prevIter();
    Urel.correctBoundaryConditions();
