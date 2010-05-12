bool closedVolume = pd.needReference();

rho = thermo->rho();

volScalarField rUA = 1.0/UEqn.A();
U = rUA*UEqn.H();

phi =
    fvc::interpolate(rho)
   *(
        (fvc::interpolate(U) & mesh.Sf())
      + fvc::ddtPhiCorr(rUA, rho, U, phi)
    )
  - fvc::interpolate(rho*rUA*gh)*fvc::snGrad(rho)*mesh.magSf();

for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
{
    fvScalarMatrix pdEqn
    (
        fvm::ddt(psi, pd)
      + fvc::ddt(psi)*pRef
      + fvc::ddt(psi, rho)*gh
      + fvc::div(phi)
      - fvm::laplacian(rho*rUA, pd)
    );

    if (corr == nCorr-1 && nonOrth == nNonOrthCorr)
    {
        pdEqn.solve(mesh.solver(pd.name() + "Final"));
    }
    else
    {
        pdEqn.solve(mesh.solver(pd.name()));
    }

    if (nonOrth == nNonOrthCorr)
    {
        phi += pdEqn.flux();
    }
}

p == pd + rho*gh + pRef;
dpdt = fvc::ddt(p);

#include "rhoEqn.H"
#include "compressibleContinuityErrs.H"

U -= rUA*(fvc::grad(pd) + fvc::grad(rho)*gh);
U.correctBoundaryConditions();


// For closed-volume cases adjust the pressure and density levels
// to obey overall mass continuity
if (closedVolume)
{
    p += (initialMass - fvc::domainIntegrate(thermo->psi()*p))
        /fvc::domainIntegrate(thermo->psi());
    rho = thermo->rho();
}
