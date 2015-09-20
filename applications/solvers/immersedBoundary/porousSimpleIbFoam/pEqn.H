    p.boundaryField().updateCoeffs();

    volScalarField AU = UEqn().A();

    U = UEqn().H()/AU;
    // Immersed boundary update
    U.correctBoundaryConditions();

    UEqn.clear();
    phi = faceIbMask*(fvc::interpolate(U) & mesh.Sf());

    // Adjust immersed boundary fluxes
    immersedBoundaryAdjustPhi(phi, U);
    adjustPhi(phi, U, p);

    // Non-orthogonal pressure corrector loop
    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(1.0/AU, p) == fvc::div(phi)
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

#   include "immersedBoundaryContinuityErrs.H"

    // Explicitly relax pressure for momentum corrector
    p.relax();

    // Momentum corrector
    U -= fvc::grad(p)/AU;
    U.correctBoundaryConditions();

