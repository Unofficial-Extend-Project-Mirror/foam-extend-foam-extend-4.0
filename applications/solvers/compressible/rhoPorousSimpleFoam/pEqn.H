if (pressureImplicitPorosity)
{
    U = trTU()&UEqn().H();
}
else
{
    U = trAU()*UEqn().H();
}

UEqn.clear();

phi = fvc::interpolate(rho*U) & mesh.Sf();
bool closedVolume = adjustPhi(phi, U, p);

while (simple.correctNonOrthogonal())
{
    tmp<fvScalarMatrix> tpEqn;

    if (pressureImplicitPorosity)
    {
        tpEqn = (fvm::laplacian(rho*trTU(), p) == fvc::div(phi));
    }
    else
    {
        tpEqn = (fvm::laplacian(rho*trAU(), p) == fvc::div(phi));
    }

    tpEqn().setReference(pRefCell, pRefValue);

    tpEqn().solve();

    if (simple.finalNonOrthogonalIter())
    {
        phi -= tpEqn().flux();
    }
}

#include "incompressible/continuityErrs.H"

// Explicitly relax pressure for momentum corrector
p.relax();

if (pressureImplicitPorosity)
{
    U -= trTU()&fvc::grad(p);
}
else
{
    U -= trAU()*fvc::grad(p);
}

U.correctBoundaryConditions();

bound(p, pMin);

// For closed-volume cases adjust the pressure and density levels
// to obey overall mass continuity
if (closedVolume)
{
    p += (initialMass - fvc::domainIntegrate(psi*p))
        /fvc::domainIntegrate(psi);
}

rho = thermo.rho();
rho.relax();
Info<< "rho max/min : " << max(rho).value() << " " << min(rho).value() << endl;
