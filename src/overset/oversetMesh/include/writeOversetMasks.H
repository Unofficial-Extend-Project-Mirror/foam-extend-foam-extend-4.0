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
        dimless,
        "zeroGradient"
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

    // Update boundary values, making sure that we skip the overset patch
    volScalarField::GeometricBoundaryField& regionIndexb =
        regionIndex.boundaryField();

    forAll(regionIndexb, patchI)
    {
        // Get the patch field
        fvPatchScalarField& ripf = regionIndexb[patchI];

        if (!isA<oversetFvPatchScalarField>(ripf))
        {
            ripf = ripf.patchInternalField();
        }
    }

    // Write region index
    regionIndex.write();
}
