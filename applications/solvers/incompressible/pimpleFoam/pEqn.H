U = rUA*UEqn().H();

if (pimple.nCorrPISO() <= 1)
{
    UEqn.clear();
}

phi = (fvc::interpolate(U) & mesh.Sf())
    + fvc::ddtPhiCorr(rUA, U, phi);

adjustPhi(phi, U, p);

// Non-orthogonal pressure corrector loop
while (pimple.correctNonOrthogonal())
{
    // Pressure corrector
    fvScalarMatrix pEqn
    (
        fvm::laplacian(rUA, p) == fvc::div(phi)
    );

    pEqn.setReference(pRefCell, pRefValue);
    pEqn.solve
    (
        mesh.solutionDict().solver(p.select(pimple.finalInnerIter()))
    );

    if (pimple.finalNonOrthogonalIter())
    {
        phi -= pEqn.flux();
    }
}

#include "continuityErrs.H"

// Explicitly relax pressure for momentum corrector except for last corrector
if (!pimple.finalIter())
{
    p.relax();
}

U -= rUA*fvc::grad(p);
U.correctBoundaryConditions();
