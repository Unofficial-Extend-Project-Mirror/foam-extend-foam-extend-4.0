// Convection-diffusion matrx
fvVectorMatrix HUEqn
(
    fvm::div(phi, U)
  + turbulence->divDevReff(U)
);

// Time derivative matrix
fvVectorMatrix ddtUEqn(fvm::ddt(U));

// Get under-relaxation factor
scalar UUrf = mesh.solutionDict().relaxationFactor(U.name());

if (oCorr == nOuterCorr - 1)
{
    if (mesh.solutionDict().relax("UFinal"))
    {
        UUrf = mesh.solutionDict().relaxationFactor("UFinal");
    }
    else
    {
        UUrf = 1;
    }
}

if (momentumPredictor)
{
    solve
    (
        ddtUEqn
      + relax(HUEqn, UUrf)
     ==
      - fvc::grad(p)
    );
}
else
{
    U = (ddtUEqn.H() + HUEqn.H() - fvc::grad(p))/(HUEqn.A() + ddtUEqn.A());
    U.correctBoundaryConditions();
}
