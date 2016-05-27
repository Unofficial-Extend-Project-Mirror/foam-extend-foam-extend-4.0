{
    volScalarField kappaEff
    (
        "kappaEff",
        turbulence->nu()/Pr + turbulence->nut()/Prt
    );

    fvScalarMatrix TEqn
    (
        fvm::div(phi, T)
      - fvm::Sp(fvc::div(phi), T)
      - fvm::laplacian(kappaEff, T)
    );

    TEqn.relax();

    TEqn.solve();

    rhok = 1.0 - beta*(T - TRef);
}
