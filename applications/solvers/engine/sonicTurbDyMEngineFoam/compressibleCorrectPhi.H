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

    pcorr == p;
    pcorr.oldTime() == p.oldTime();

    phi = fvc::interpolate(rho)
        *((fvc::interpolate(U) & mesh.Sf()) - fvc::meshPhi(rho, U)); 

    for(int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pcorrEqn
        (
            fvm::ddt(psi, pcorr)
          + fvc::div(phi)
          - fvm::laplacian(rho*rUA, pcorr)
        );

        pcorrEqn.solve();

        if (nonOrth == nNonOrthCorr)
        {
            phi += pcorrEqn.flux(); 
        }
    }
}
