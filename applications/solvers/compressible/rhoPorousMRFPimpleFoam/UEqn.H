// Solve the Momentum equation

tmp<fvVectorMatrix> UEqn
(
    pZones.ddt(rho, U)
  + fvm::div(phi, U)
  + turbulence->divDevRhoReff()
);

UEqn().relax
(
    mesh.solutionDict().equationRelaxationFactor(U.select(pimple.finalIter()))
);

mrfZones.addCoriolis(rho, UEqn());
pZones.addResistance(UEqn());

volScalarField rUA = 1.0/UEqn().A();

if (pimple.momentumPredictor())
{
    solve
    (
        UEqn() == -fvc::grad(p),
        mesh.solutionDict().solver((U.select(pimple.finalIter())))
    );
}
else
{
    U = rUA*(UEqn().H() - fvc::grad(p));
    U.correctBoundaryConditions();
}
