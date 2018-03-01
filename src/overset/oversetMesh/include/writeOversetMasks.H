{
    const oversetMesh& om = oversetMesh::New(mesh);

    om.gamma().write();
    om.gammaExt().write();
    om.sGamma().write();
    om.oversetTypes().write();

    // Create region ID field
    volScalarField regionIndex
    (
        IOobject
        (
            "regionIndex",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimless
    );
    scalarField& regionIndexIn = regionIndex.internalField();

    forAll (om.regions(), regionI)
    {
        const cellZone& cz = om.regions()[regionI].zone();

        forAll (cz, cellI)
        {
            regionIndexIn[cz[cellI]] = regionI;
        }
    }

    // Update boundary values
    forAll (regionIndex.boundaryField(), patchI)
    {
        if (!regionIndex.boundaryField()[patchI].coupled())
        {
            regionIndex.boundaryField()[patchI] =
                regionIndex.boundaryField()[patchI].patchInternalField();
        }
    }

//     regionIndex.write();
}
