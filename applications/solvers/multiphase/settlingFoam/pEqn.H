volScalarField rUA = 1.0/UEqn.A();

surfaceScalarField rUAf
(
    "(rho*(1|A(U)))",
    fvc::interpolate(rho)*fvc::interpolate(rUA)
);

U = rUA*UEqn.H();
phi =
    fvc::interpolate(rho)
   *(
       (fvc::interpolate(U) & mesh.Sf())
     + fvc::ddtPhiCorr(rUA, rho, U, phi)
    );

surfaceScalarField phiU("phiU", phi);
phi -= ghf*fvc::snGrad(rho)*rUAf*mesh.magSf();

while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix pEqn
    (
        fvm::laplacian(rUAf, p) == fvc::ddt(rho) + fvc::div(phi)
    );

    pEqn.setReference(pRefCell, pRefValue);
    pEqn.solve();

    if (pimple.finalNonOrthogonalIter())
    {
        phi -= pEqn.flux();
    }
}

#include "rhoEqn.H"
#include "compressibleContinuityErrs.H"

U += rUA*fvc::reconstruct((phi - phiU)/rUAf);
U.correctBoundaryConditions();
