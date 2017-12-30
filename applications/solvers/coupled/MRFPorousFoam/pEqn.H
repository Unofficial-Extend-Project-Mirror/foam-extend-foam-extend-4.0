// Pressure parts of the continuity equation
tmp<fvScalarMatrix> tpEqn;
tmp<surfaceScalarField> tpresSource;

if (addMRF || addPorosity)
{
    // Collect pressure source with tensorial 1/Ap
    const volTensorField& rTU = trTU();

    tpresSource =
    (
        (mesh.Sf() & fvc::interpolate(rTU))
      & fvc::interpolate(fvc::grad(p))
    );
    const surfaceScalarField& presSource = tpresSource();

    // Assemble pressure matrix with tensorial 1/Ap
    tpEqn =
    (
      - fvm::laplacian(rTU, p)
     ==
      - fvc::div(presSource)
    );
}
else
{
    // Collect pressure source with scalar 1/Ap
    const volScalarField& rAU = trAU();

    tpresSource =
    (
        fvc::interpolate(rAU)*
        (mesh.Sf() & fvc::interpolate(fvc::grad(p, "grad(pSource)")))
    );
    const surfaceScalarField& presSource = tpresSource();

    // Assemble pressure matrix with tensorial 1/Ap
    tpEqn =
    (
        - fvm::laplacian(rAU, p)
       ==
        - fvc::div(presSource)
    );
}

tpEqn().setReference(pRefCell, pRefValue);

UpEqn.insertEquation(3, tpEqn());
