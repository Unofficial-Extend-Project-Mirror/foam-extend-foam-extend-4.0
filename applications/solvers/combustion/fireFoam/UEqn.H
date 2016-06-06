fvVectorMatrix UEqn
(
    fvm::ddt(rho, U)
  + fvm::div(phi, U)
  + turbulence->divDevRhoReff()
);

UEqn.relax();

solve
(
    UEqn
  ==
    fvc::reconstruct
    (
        fvc::interpolate(rho)*(g & mesh.Sf())
      - fvc::snGrad(p)*mesh.magSf()
    ),
    mesh.solutionDict().solver((U.select(pimple.finalIter())))
);
