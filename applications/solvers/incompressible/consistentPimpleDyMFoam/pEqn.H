{
    p.boundaryField().updateCoeffs();

    // Prepare clean Ap without time derivative contribution and
    // without contribution from under-relaxation
    // HJ, 26/Oct/2015
    aU = HUEqn.A();

    // Store velocity under-relaxation point before using U for
    // the flux precursor
    U.storePrevIter();

    U = HUEqn.H()/aU;
    phi = (fvc::interpolate(U) & mesh.Sf());

    // Global flux balance
    adjustPhi(phi, U, p);

    for (int nonOrth = 0; nonOrth <= nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(1/aU, p) == fvc::div(phi)
        );

        pEqn.setReference(pRefCell, pRefValue);

        if
        (
//             oCorr == nOuterCorr - 1
            corr == nCorr - 1
         && nonOrth == nNonOrthCorr
        )
        {
            pEqn.solve
            (
                mesh.solutionDict().solver(p.name() + "Final")
            );
        }
        else
        {
            pEqn.solve(mesh.solutionDict().solver(p.name()));
        }

        if (nonOrth == nNonOrthCorr)
        {
            phi -= pEqn.flux();
        }
    }

    // Explicitly relax pressure for momentum corrector
    if (oCorr != nOuterCorr - 1)
    {
        p.relax();
    }

    // Make the fluxes relative to the mesh motion
    fvc::makeRelative(phi, U);

#   include "movingMeshContinuityErrs.H"

    U = UUrf*
        (
            1.0/(aU + ddtUEqn.A())*
            (
                U*aU - fvc::grad(p) + ddtUEqn.H()
            )
        )
      + (1 - UUrf)*U.prevIter();
    U.correctBoundaryConditions();
}
