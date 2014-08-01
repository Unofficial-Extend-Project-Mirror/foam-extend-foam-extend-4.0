// Momentum equation
fvVectorMatrix UEqn
(
    fvm::div(phi, U)
  + fvc::div(-phi, U, "div(phi,U)")
  + turbulence->divDevReff(U)
);

UEqn.relax();

UpEqn.insertEquation(0, UEqn);
UpEqn.insertPicardTensor(0, U, phi);
