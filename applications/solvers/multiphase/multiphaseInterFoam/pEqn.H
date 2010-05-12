{
    volScalarField rUA = 1.0/UEqn.A();
    surfaceScalarField rUAf = fvc::interpolate(rUA);

    U = rUA*UEqn.H();

    surfaceScalarField phiU
    (
        "phiU",
        (fvc::interpolate(U) & mesh.Sf()) + fvc::ddtPhiCorr(rUA, rho, U, phi)
    );

    phi = phiU +
        (
            mixture.surfaceTensionForce()
          - ghf*fvc::snGrad(rho)
        )*rUAf*mesh.magSf();

    adjustPhi(phi, U, pd);

    for(int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pdEqn
        (
            fvm::laplacian(rUAf, pd) == fvc::div(phi)
        );

        pdEqn.setReference(pdRefCell, pdRefValue);

        if (corr == nCorr-1)
        {
            pdEqn.solve(mesh.solver(pd.name() + "Final"));
        }
        else
        {
            pdEqn.solve(mesh.solver(pd.name()));
        }

        if (nonOrth == nNonOrthCorr)
        {
            phi -= pdEqn.flux();
        }
    }

    U += rUA*fvc::reconstruct((phi - phiU)/rUAf);
    U.correctBoundaryConditions();
}
