    // Pressure parts of the continuity equation
    surfaceScalarField presSource
    (
         "presSource",
         fvc::interpolate(rAU)*
         (fvc::interpolate(fvc::grad(p)) & mesh.Sf())
    );

    fvScalarMatrix pEqn
    (
      - fvm::laplacian(rAU, p)
     ==
      - fvc::div(presSource)
    );

    pEqn.setReference(pRefCell, pRefValue);

    UpEqn.insertEquation(3, pEqn);
