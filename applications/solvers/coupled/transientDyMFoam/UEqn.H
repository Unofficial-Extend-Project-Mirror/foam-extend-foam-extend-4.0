{
    // Momentum equation
    fvVectorMatrix UEqn
    (
        fvm::ddt(U)
      + fvm::div(phi, U)
      + turbulence->divDevReff()
    );

    rAU = 1.0/UEqn.A();

    // Under-relax momentum. Note this will destroy the H and A
    UEqn.relax();

    // Insert momentum equation
    UpEqn.insertEquation(0, UEqn);

#   include "addBlockCoupledBC.H"
}
