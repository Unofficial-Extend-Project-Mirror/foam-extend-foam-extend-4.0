{
    volScalarField rUA = 1.0/UEqn.A();
    surfaceScalarField rUAf = fvc::interpolate(rUA);

    U = rUA*UEqn.H();

    // Update boundary velocity for consistency with the flux
    mrfZones.correctBoundaryVelocity(U);

    surfaceScalarField phiU
    (
        "phiU",
        (fvc::interpolate(U) & mesh.Sf())
    //+ fvc::ddtPhiCorr(rUA, rho, U, phi)
    );
    mrfZones.relativeFlux(phiU);

    phi = phiU +
        (
            fvc::interpolate(interface.sigmaK())*
            fvc::snGrad(alpha1)*mesh.magSf()
          + fvc::interpolate(rho)*(g & mesh.Sf())
        )*rUAf;
    adjustPhi(phi, U, p);

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pdEqn
        (
            fvm::laplacian(rUAf, pd) == fvc::div(phi)
        );

        pdEqn.setReference(pdRefCell, pdRefValue);

        pdEqn.solve
        (
            mesh.solutionDict().solver(pd.select(pimple.finalInnerIter()))
        );

        if (pimple.finalNonOrthogonalIter())
        {
            phi -= pdEqn.flux();
        }
    }

    U += rUA*fvc::reconstruct((phi - phiU)/rUAf);
    U.correctBoundaryConditions();
}
