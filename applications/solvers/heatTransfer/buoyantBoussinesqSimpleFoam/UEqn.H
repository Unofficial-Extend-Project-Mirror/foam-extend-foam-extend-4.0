    // Solve the momentum equation

    tmp<fvVectorMatrix> UEqn
    (
        fvm::div(phi, U)
      + turbulence->divDevReff()
    );

    solve
    (
        relax(UEqn())
      ==
       -fvc::grad(p)
      + fvc::reconstruct(fvc::interpolate(rhok)*(g & mesh.Sf()))
    );
