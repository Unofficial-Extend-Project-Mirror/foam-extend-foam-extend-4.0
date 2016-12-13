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

    // Insert the additional components. Note this will destroy the H and A

    UEqn += fvm::SuSp(-divPhi, U) + divPhi*U;
    UEqn.relax();

    UpEqn.insertEquation(0, UEqn);
}
