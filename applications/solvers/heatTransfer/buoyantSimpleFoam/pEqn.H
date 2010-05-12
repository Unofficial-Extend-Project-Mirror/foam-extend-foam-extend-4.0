volScalarField rUA = 1.0/UEqn().A();
U = rUA*UEqn().H();
UEqn.clear();
phi = fvc::interpolate(rho)*(fvc::interpolate(U) & mesh.Sf());
bool closedVolume = adjustPhi(phi, U, p);
phi -= fvc::interpolate(rho*gh*rUA)*fvc::snGrad(rho)*mesh.magSf();

for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
{
    fvScalarMatrix pdEqn
    (
        fvm::laplacian(rho*rUA, pd) == fvc::div(phi)
    );

    pdEqn.setReference(pdRefCell, pdRefValue);
    // retain the residual from the first iteration
    if (nonOrth == 0)
    {
        eqnResidual = pdEqn.solve().initialResidual();
        maxResidual = max(eqnResidual, maxResidual);
    }
    else
    {
        pdEqn.solve();
    }

    if (nonOrth == nNonOrthCorr)
    {
        phi -= pdEqn.flux();
    }
}

#include "continuityErrs.H"

// Explicitly relax pressure for momentum corrector
pd.relax();

p = pd + rho*gh + pRef;

U -= rUA*(fvc::grad(pd) + fvc::grad(rho)*gh);
U.correctBoundaryConditions();

// For closed-volume cases adjust the pressure and density levels
// to obey overall mass continuity
if (closedVolume)
{
    p += (initialMass - fvc::domainIntegrate(thermo->psi()*p))
        /fvc::domainIntegrate(thermo->psi());
}

rho = thermo->rho();
rho.relax();
Info<< "rho max/min : " << max(rho).value() << " " << min(rho).value() << endl;
