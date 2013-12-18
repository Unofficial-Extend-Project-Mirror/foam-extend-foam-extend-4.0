    porousZones pZones(mesh);

    // Creates the porosity field for MULES
    volScalarField porosity
    (
        IOobject
        (
            "porosity",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("NULL", dimless, 1.0),
        "zeroGradient"
    );

    forAll( pZones, zoneI )
    {
        const label & zoneId( pZones[zoneI].zoneId() );

        const labelList & cells(mesh.cellZones()[zoneId]);

        const scalar & zonePorosity( pZones[zoneI].porosity() );

        forAll( cells, cellI )
        {
            porosity[cells[cellI]] = zonePorosity;
        }
    }

    porosity.write();
