volScalarField AU = UEqn().A();
U = UEqn().H()/AU;
UEqn.clear();
phi = fvc::interpolate(rho*U) & mesh.Sf();
bool closedVolume = adjustPhi(phi, U, p);

for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
{
    fvScalarMatrix pEqn
    (
        fvm::laplacian(rho/AU, p) == fvc::div(phi)
    );

    pEqn.setReference(pRefCell, pRefValue);
    // retain the residual from the first iteration
    if (nonOrth == 0)
    {
        eqnResidual = pEqn.solve().initialResidual();
        maxResidual = max(eqnResidual, maxResidual);
    }
    else
    {
        pEqn.solve();
    }

    if (nonOrth == nNonOrthCorr)
    {
        phi -= pEqn.flux();
    }
}

#include "incompressible/continuityErrs.H"

// Explicitly relax pressure for momentum corrector
p.relax();

U -= fvc::grad(p)/AU;
U.correctBoundaryConditions();

bound(p, pMin);

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
