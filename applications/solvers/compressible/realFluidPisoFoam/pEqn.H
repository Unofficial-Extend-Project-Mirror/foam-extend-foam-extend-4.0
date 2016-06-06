rho = thermo.rho();

volScalarField rUA = 1.0/UEqn.A();
U = rUA*UEqn.H();

phi =
    fvc::interpolate(rho)*
    (
        (fvc::interpolate(U) & mesh.Sf())
      + fvc::ddtPhiCorr(rUA, rho, U, phi)
    );

while (piso.correctNonOrthogonal())
{
    fvScalarMatrix pEqn
    (
        psi*fvm::ddt(p)
      + drhodh*fvc::ddt(h)
      + fvc::div(phi)
      - fvm::laplacian(rho*rUA, p)
    );

    pEqn.solve();

    if (piso.finalNonOrthogonalIter())
    {
        phi += pEqn.flux();
    }
}

#include "rhoEqn.H"
#include "compressibleContinuityErrs.H"

U -= rUA*fvc::grad(p);
U.correctBoundaryConditions();

DpDt = fvc::DDt(surfaceScalarField("phiU", phi/fvc::interpolate(rho)), p);
