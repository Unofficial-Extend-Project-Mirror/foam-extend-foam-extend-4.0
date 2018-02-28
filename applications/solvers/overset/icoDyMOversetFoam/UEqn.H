// Momentum equation
fvVectorMatrix UEqn
(
    fvm::ddt(U)
  + fvm::div(phi, U)
  - fvm::laplacian(nu, U)
);

if (piso.momentumPredictor())
{
    solve(UEqn == -fvc::grad(p));
}
