// Momentum equation
fvVectorMatrix UEqn
(
    fvm::div(phi, U)
  + turbulence->divDevReff(U)
);

UEqn.relax();

UpEqn.insertEquation(0, UEqn);
