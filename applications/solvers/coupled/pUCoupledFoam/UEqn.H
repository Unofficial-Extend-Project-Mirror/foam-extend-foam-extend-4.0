{
    volScalarField divPhi
    (
        "divPhi",
        fvc::div(phi)
    );

    // Momentum equation
    fvVectorMatrix UEqn
    (
        fvm::div(phi, U)
      + turbulence->divDevReff()
    );

    rAU = 1.0/UEqn.A();

    // Under-relax momentum. Note this will destroy the H and A
    UEqn.relax();

    UpEqn.insertEquation(0, UEqn);
}
