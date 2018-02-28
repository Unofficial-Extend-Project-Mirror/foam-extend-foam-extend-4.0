{
    wordList pcorrTypes(p.boundaryField().types());

    for (label i = 0; i < p.boundaryField().size(); i++)
    {
        if (p.boundaryField()[i].fixesValue())
        {
            pcorrTypes[i] = fixedValueFvPatchScalarField::typeName;
        }
    }

    volScalarField pcorr
    (
        IOobject
        (
            "pcorr",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("pcorr", p.dimensions(), 0),
        pcorrTypes
    );

    // Initialise flux with interpolated velocity
    phi = faceOversetMask*(fvc::interpolate(U) & mesh.Sf());

    // Adjust overset fluxes
    oversetAdjustPhi(phi, U); // Fringe flux adjustment
    globalOversetAdjustPhi(phi, U, pcorr); // Global flux adjustment

    while (piso.correctNonOrthogonal())
    {
        fvScalarMatrix pcorrEqn
        (
            fvm::laplacian(rAU, pcorr) == fvc::div(phi)
        );

        // Adjust non-orthogonal fluxes if necessary
        om.correctNonOrthoFluxes(pcorrEqn, U);

        pcorrEqn.setReference(pRefCell, pRefValue);
        pcorrEqn.solve();

        if (piso.finalNonOrthogonalIter())
        {
            phi -= pcorrEqn.flux();
        }

        // Fluxes are corrected to absolute velocity and further corrected
        // later.  HJ, 6/Feb/2009
    }

#   include "oversetContinuityErrs.H"
}
