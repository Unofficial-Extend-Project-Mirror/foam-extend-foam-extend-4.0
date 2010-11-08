{
    surfaceScalarField rUAf = 1.0/fvc::interpolate(UEqn.A());

    U = UEqn.H()/UEqn.A();

    phi = fvc::interpolate(U) & mesh.Sf();

    // Bug fix: must change name of phi on copy to keep objectRegistry happy
    // HJ, 7/Nov/2010
    surfaceScalarField phiU
    (
        "phiU",
        phi
    );

    surfaceScalarField phip = fvc::interpolate(psiByRho)*phi;

    for(int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(rUAf, p)
        );

        solve
        (
             psiByRho*fvm::ddt(p)
           + fvm::div(phip, p) - fvm::Sp(fvc::div(phip), p)
           + fvc::div(phi) - pEqn
        );

        if (nonOrth == nNonOrthCorr)
        {
            phi -= pEqn.flux();
        }
    }

    U += fvc::reconstruct(phi - phiU);
    U.correctBoundaryConditions();
}
