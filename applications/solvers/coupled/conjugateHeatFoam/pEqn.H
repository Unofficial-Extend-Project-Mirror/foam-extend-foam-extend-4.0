{
    volScalarField rUA("rUA", 1.0/UEqn.A());
    surfaceScalarField rUAf("(1|A(U))", fvc::interpolate(rUA));

    U = rUA*UEqn.H();

    phi = (fvc::interpolate(U) & mesh.Sf())
        + fvc::ddtPhiCorr(rUA, U, phi);

    surfaceScalarField buoyancyPhi(rUAf*ghf*fvc::snGrad(rhok)*mesh.magSf());
    phi -= buoyancyPhi;

    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix p_rghEqn
        (
            fvm::laplacian(rUAf, p_rgh) == fvc::div(phi)
        );

        p_rghEqn.setReference(pRefCell, pRefValue);

        if (corr == nCorr-1 && nonOrth == nNonOrthCorr)
        {
            p_rghEqn.solve(mesh.solutionDict().solver(p_rgh.name() + "Final"));
        }
        else
        {
            p_rghEqn.solve(mesh.solutionDict().solver(p_rgh.name()));
        }

        if (nonOrth == nNonOrthCorr)
        {
            // Calculate the conservative fluxes
            phi -= p_rghEqn.flux();

            // Explicitly relax pressure for momentum corrector
            p_rgh.relax();

            // Correct the momentum source with the pressure gradient flux
            // calculated from the relaxed pressure
            U -= rUA*fvc::reconstruct((buoyancyPhi + p_rghEqn.flux())/rUAf);
            U.correctBoundaryConditions();
        }
    }

    #include "continuityErrs.H"
}
