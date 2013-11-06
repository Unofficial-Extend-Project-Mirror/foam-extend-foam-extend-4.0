    Info << "Reading coupling properties" << endl;
    IOdictionary couplingProperties
    (
        IOobject
        (
            "couplingProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // Read solid patch data
    word solidPatchName(couplingProperties.lookup("solidPatch"));

    label solidPatchID =
        stressMesh.boundaryMesh().findPatchID(solidPatchName);


    // Read fluid patch data
    word fluidPatchName(couplingProperties.lookup("fluidPatch"));

    label fluidPatchID =
        mesh.boundaryMesh().findPatchID(fluidPatchName);


    if (solidPatchID < 0 || fluidPatchID < 0)
    {
        FatalErrorIn(args.executable())
            << "Problem with patch interpolation definition"
            << abort(FatalError);
    }


    // Create interpolators
    patchToPatchInterpolation interpolatorFluidSolid
    (
        mesh.boundaryMesh()[fluidPatchID],
        stressMesh.boundaryMesh()[solidPatchID]
    );

    patchToPatchInterpolation interpolatorSolidFluid
    (
        stressMesh.boundaryMesh()[solidPatchID],
        mesh.boundaryMesh()[fluidPatchID]
    );


    // Grab solid patch field
    tractionDisplacementFvPatchVectorField& tForce =
        refCast<tractionDisplacementFvPatchVectorField>
        (
            Usolid.boundaryField()[solidPatchID]
        );


    // Grab motion field

    // Read fluid patch data
    word movingRegionName(couplingProperties.lookup("movingRegion"));

    const fvMesh& motionMesh =
        runTime.objectRegistry::lookupObject<fvMesh>(movingRegionName);

    tetPointVectorField& motionU =
        const_cast<tetPointVectorField&>
        (
            motionMesh.objectRegistry::lookupObject<tetPointVectorField>
            (
                "motionU"
            )
        );


    fixedValueTetPolyPatchVectorField& motionUFluidPatch =
        refCast<fixedValueTetPolyPatchVectorField>
        (
            motionU.boundaryField()[fluidPatchID]
        );

    tetPolyPatchInterpolation tppi
    (
        refCast<const faceTetPolyPatch>(motionUFluidPatch.patch())
    );
