{
    p.boundaryField().updateCoeffs();

    rAU = 1.0/UEqn.A();
    oversetFvPatchScalarField::oversetInterpolate(rAU); // Overset update
    surfaceScalarField rAUf = fvc::interpolate(rAU);

    U = rAU*UEqn.H();
    oversetFvPatchVectorField::oversetInterpolate(U); // Overset update

    surfaceScalarField phiU
    (
        "phiU",
        fvc::interpolate(U) & mesh.Sf()
    );

    phi = phiU;

    // Adjust overset fluxes
    oversetAdjustPhi(phi, U); // Fringe flux adjustment
    globalOversetAdjustPhi(phi, U, p); // Global flux adjustment

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(rAUf, p) == fvc::div(phi)
        );

        // Adjust non-orthogonal fringe fluxes if necessary
        om.correctNonOrthoFluxes(pEqn, U);

        pEqn.setReference(pRefCell, pRefValue);

        if
        (
            pimple.corrPISO() == pimple.nCorrPISO()
         && pimple.finalNonOrthogonalIter()
        )
        {
            pEqn.solve(mesh.solutionDict().solver("pFinal"));
        }
        else
        {
            pEqn.solve();
        }

        if (pimple.finalNonOrthogonalIter())
        {
            phi -= pEqn.flux();
        }

        // Perform overset interpolation (after flux reconstruction)
        oversetFvPatchScalarField::oversetInterpolate(p);
    }

#   include "oversetContinuityErrs.H"

    // Explicitly relax pressure for momentum corrector
    if (!pimple.finalIter())
    {
        p.relax();
    }

#   include "movingMeshContinuityErrs.H"

    U += rAU*fvc::reconstruct((phi - phiU)/rAUf);
    U.correctBoundaryConditions();
    oversetFvPatchVectorField::oversetInterpolate(U); // Overset update
    // Note: if implicit fringe is switched on, we are doing the interpolation
    // twice (once in correctBoundaryConditions and once in oversetInterpolate)
    // Reorganize. VV, 4/Oct/2016.

    // Make the fluxes relative to the mesh motion
    fvc::makeRelative(phi, U);
}
