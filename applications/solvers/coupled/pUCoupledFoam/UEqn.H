// Momentum equation
tmp<fvVectorMatrix> UEqn
(
    fvm::ddt(U)
  + fvm::div(phi, U)
  + turbulence->divDevReff(U)
);

UEqn().relax();

UpEqn.insertEquation(0, UEqn());
