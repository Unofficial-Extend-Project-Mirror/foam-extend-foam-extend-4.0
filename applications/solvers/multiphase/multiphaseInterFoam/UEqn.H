    surfaceScalarField muf = mixture.muf();

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(mixture.rhoPhi(), U)
      - fvm::laplacian(muf, U)
      - (fvc::grad(U) & fvc::grad(muf))
    //- fvc::div(muf*(mesh.Sf() & fvc::interpolate(T(fvc::grad(U)))))
    );
