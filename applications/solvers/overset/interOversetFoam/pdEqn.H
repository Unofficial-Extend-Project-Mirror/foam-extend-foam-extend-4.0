{
    pd.boundaryField().updateCoeffs();

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

    phi = phiU
      + (
            fvc::interpolate(interface.sigmaK())*fvc::snGrad(alpha1)
          - ghf*fvc::snGrad(rho)
        )*rAUf*mesh.magSf();

    // Adjust overset fluxes
    oversetAdjustPhi(phi, U); // Fringe flux adjustment
    globalOversetAdjustPhi(phi, U, p); // Global flux adjustment

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pdEqn
        (
            fvm::laplacian(rAUf, pd) == fvc::div(phi)
        );

        // Adjust non-orthogonal fringe fluxes if necessary
        om.correctNonOrthoFluxes(pdEqn, U);

        pdEqn.setReference(pRefCell, pRefValue);

        pdEqn.solve
        (
            mesh.solutionDict().solver(pd.select(pimple.finalInnerIter()))
        );

        if (pimple.finalNonOrthogonalIter())
        {
            phi -= pdEqn.flux();
        }

        // Perform overset interpolation (after flux reconstruction)
        oversetFvPatchScalarField::oversetInterpolate(p);
    }

    // Explicitly relax pressure for momentum corrector
    if (!pimple.finalIter())
    {
        p.relax();
    }

    U += rAU*fvc::reconstruct((phi - phiU)/rAUf);
    U.correctBoundaryConditions();
    oversetFvPatchVectorField::oversetInterpolate(U); // Overset update
    // Note: if implicit fringe is switched on, we are doing the interpolation
    // twice (once in correctBoundaryConditions and once in oversetInterpolate)
    // Reorganize. VV, 4/Oct/2016.
}
