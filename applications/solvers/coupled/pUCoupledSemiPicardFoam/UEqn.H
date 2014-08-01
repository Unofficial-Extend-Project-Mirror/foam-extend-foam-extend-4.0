// Momentum equation
fvVectorMatrix UEqn
(
    fvm::div(2*phi, U,"div(phi,U)")
  + fvc::div(-phi, U,"div(phi,U)")
  + turbulence->divDevReff(U)
);

UEqn.relax();

UpEqn.insertEquation(0, UEqn);
