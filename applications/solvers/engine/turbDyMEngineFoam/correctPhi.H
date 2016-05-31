{
    wordList pcorrTypes(p.boundaryField().types());

    for (label i=0; i<p.boundaryField().size(); i++)
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
        dimensionedScalar("pcorr", p.dimensions(), 0.0),
        pcorrTypes
    );

    phi = fvc::interpolate(U) & mesh.Sf();

    mesh.schemesDict().setFluxRequired(pcorr.name());

    while (piso.correctNonOrthogonal())
    {
        fvScalarMatrix pcorrEqn
        (
            fvm::laplacian(rUA, pcorr) == fvc::div(phi)
        );

        pcorrEqn.setReference(pRefCell, pRefValue);
        pcorrEqn.solve();

        if (piso.finalNonOrthogonalIter())
        {
            phi -= pcorrEqn.flux();
        }
    }
}

#include "continuityErrs.H"
