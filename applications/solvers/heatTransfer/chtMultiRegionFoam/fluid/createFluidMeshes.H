    PtrList<fvMesh> fluidRegions(rp.fluidRegionNames().size());

    forAll(rp.fluidRegionNames(), i)
    {
        Info<< "Create fluid mesh for region " << rp.fluidRegionNames()[i]
            << " for time = " << runTime.timeName() << nl << endl;

        fluidRegions.set
        (
            i,
            new fvMesh
            (
                IOobject
                (
                    rp.fluidRegionNames()[i],
                    runTime.timeName(),
                    runTime,
                    IOobject::MUST_READ
                )
            )
        );

        // Force calculation of geometric properties to prevent it being done
        // later in e.g. some boundary evaluation
        //(void)fluidRegions[i].weights();
        //(void)fluidRegions[i].deltaCoeffs();

        // Attach environmental properties to each region
        autoPtr<IOdictionary> environmentalProperties
        (
            new IOdictionary
            (
                IOobject
                (
                    "environmentalProperties",
                    runTime.constant(),
                    fluidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            )
        );

        environmentalProperties.ptr()->store();
    }
