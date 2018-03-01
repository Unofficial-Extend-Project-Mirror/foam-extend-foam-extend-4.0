{
    p.boundaryField().updateCoeffs();

    rAU = 1.0/UEqn().A();
    oversetFvPatchScalarField::oversetInterpolate(rAU); // Overset update

    U = rAU*UEqn().H();
    oversetFvPatchVectorField::oversetInterpolate(U); // Overset update

    UEqn.clear();
    phi = (fvc::interpolate(U) & mesh.Sf());

    // Adjust overset fluxes
    oversetAdjustPhi(phi, U); // Fringe flux adjustment
    globalOversetAdjustPhi(phi, U, p); // Global flux adjustment

    // Non-orthogonal pressure corrector loop
    while (simple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(rAU, p) == fvc::div(phi)
        );

        // Adjust non-orthogonal fringe fluxes if necessary
        om.correctNonOrthoFluxes(pEqn, U);

        pEqn.setReference(pRefCell, pRefValue);

        // Retain the residual from the first iteration
        pEqn.solve();

        if (simple.finalNonOrthogonalIter())
        {
            phi -= pEqn.flux();
        }

        // Perform overset interpolation (after flux reconstruction)
        oversetFvPatchScalarField::oversetInterpolate(p);
    }

#   include "oversetContinuityErrs.H"

    // Explicitly relax pressure for momentum corrector
    p.relax();

    // Momentum corrector
    U -= rAU*fvc::grad(p);
    U.correctBoundaryConditions();
    oversetFvPatchVectorField::oversetInterpolate(U); // Overset update
    // Note: if implicit fringe is switched on, we are doing the interpolation
    // twice (once in correctBoundaryConditions and once in oversetInterpolate)
    // Reorganize. VV, 4/Oct/2016.
}
