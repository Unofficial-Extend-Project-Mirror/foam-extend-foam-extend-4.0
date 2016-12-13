bool closedVolume = false;

rho = thermo.rho();

volScalarField rUA = 1.0/UEqn.A();
surfaceScalarField rhorUAf("(rho*(1|A(U)))", fvc::interpolate(rho*rUA));
U = rUA*UEqn.H();

surfaceScalarField phiU
(
    fvc::interpolate(rho)
   *(
        (fvc::interpolate(U) & mesh.Sf())
      + fvc::ddtPhiCorr(rUA, rho, U, phi)
    )
);

phi = phiU + rhorUAf*fvc::interpolate(rho)*(g & mesh.Sf());

while (pimple.correctNonOrthogonal())
{
       surfaceScalarField rhorUAf = fvc::interpolate(rho*rUA);

        fvScalarMatrix pEqn
        (
            fvm::ddt(psi,p)
          + fvc::div(phi)
          - fvm::laplacian(rhorUAf, p)
        );

        closedVolume = p.needReference();

        pEqn.solve
        (
            mesh.solutionDict().solver(p.select(pimple.finalInnerIter()))
        );

        if (pimple.finalNonOrthogonalIter())
        {
            phi += pEqn.flux();
        }
}

DpDt = fvc::DDt(surfaceScalarField("phiU", phi/fvc::interpolate(rho)), p);

#include "rhoEqn.H"
#include "compressibleContinuityErrs.H"

U += rUA*fvc::reconstruct((phi - phiU)/rhorUAf);
U.correctBoundaryConditions();

// For closed-volume cases adjust the pressure and density levels
// to obey overall mass continuity
if (closedVolume)
{
    p +=
        (initialMass - fvc::domainIntegrate(thermo.psi()*p))
       /fvc::domainIntegrate(thermo.psi());
    rho = thermo.rho();
}
