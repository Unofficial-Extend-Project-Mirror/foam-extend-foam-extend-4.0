{
    rho = thermo.rho();

    rUA = 1.0/UEqn.A();
    H = UEqn.H();
    U = rUA*UEqn.H();

    phi = fvc::interpolate(rho)
         *((fvc::interpolate(U) & mesh.Sf()) - fvc::meshPhi(rho, U));

    // Store pressure for under-relaxation
    p.storePrevIter();

    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pEqn
        (
            fvm::ddt(psi, p)
          + fvc::div(phi)
          - fvm::laplacian(rho*rUA, p)
        );

        pEqn.solve();

        if (nonOrth == nNonOrthCorr)
        {
            phi += pEqn.flux();
        }
    }

#   include "rhoEqn.H"
#   include "compressibleContinuityErrs.H"

    // Warning:
    // rho does not carry working boundary conditions and needs to be updated
    // strictly according to the thermodynamics package
    // HJ, 22/Aug/2007
    thermo.correct();
    rho = thermo.rho();

    DpDt = fvc::DDt(surfaceScalarField("phiU", phi/fvc::interpolate(rho)), p);

    U -= rUA*fvc::grad(p);
    U.correctBoundaryConditions();
}

