{
    volScalarField pcorr("pcorr", p);
    pcorr *= 0;

    // Initialise flux with interpolated velocity
    phi = fvc::interpolate(U) & mesh.Sf();

    adjustPhi(phi, U, pcorr);

    mesh.schemesDict().setFluxRequired(pcorr.name());

    dimensionedScalar rUAf("(1|A(U))", dimTime, 1.0);

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pcorrEqn
        (
            fvm::laplacian(rUAf, pcorr) == fvc::div(phi)
        );

        pcorrEqn.setReference(pRefCell, pRefValue);
        pcorrEqn.solve();

        if (pimple.finalNonOrthogonalIter())
        {
            phi -= pcorrEqn.flux();
        }

        // Fluxes are corrected to absolute velocity and further corrected
        // later.  HJ, 6/Feb/2009
    }

    #   include "continuityErrs.H"
}
