    Info << "\nReading coupling properties" << endl;
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

    label solidPatchID
    (
        stressMesh.boundaryMesh().findPatchID(solidPatchName)
    );

    if (solidPatchID < 0)
    {
        FatalErrorIn(args.executable())
            << "Problem with finding solid patch"
                << abort(FatalError);
    }

    word solidZoneName(couplingProperties.lookup("solidZone"));

    label solidZoneID =
        stressMesh.faceZones().findZoneID(solidZoneName);

    if (solidZoneID < 0)
    {
        FatalErrorIn(args.executable())
            << "Problem with finding solid zone"
                << abort(FatalError);
    }

//     const label solidPatchStart =
//         stressMesh.boundaryMesh()[solidPatchID].start();

//     const faceZone& solidZone = stressMesh.faceZones()[solidZoneID];


    // Read fluid patch data
    word fluidPatchName(couplingProperties.lookup("fluidPatch"));

    label fluidPatchID =
        mesh.boundaryMesh().findPatchID(fluidPatchName);

    if (fluidPatchID < 0)
    {
        FatalErrorIn(args.executable())
            << "Problem with finding fluid patch"
            << abort(FatalError);
    }

    word fluidZoneName(couplingProperties.lookup("fluidZone"));

    label fluidZoneID =
        mesh.faceZones().findZoneID(fluidZoneName);

    if (fluidZoneID < 0)
    {
        FatalErrorIn(args.executable())
            << "Problem with finding fluid zone"
                << abort(FatalError);
    }

//     const label fluidPatchStart =
//         mesh.boundaryMesh()[fluidPatchID].start();

//     const faceZone& fluidZone = mesh.faceZones()[fluidZoneID];


    // Check mesh motion solver type
    bool feMotionSolver =
        mesh.objectRegistry::foundObject<tetPointVectorField>
        (
            "motionU"
        );
    bool fvMotionSolver =
        mesh.objectRegistry::foundObject<pointVectorField>
        (
            "pointMotionU"
        );


    // Grab solid patch field
    if
    (
        DU.boundaryField()[solidPatchID].type()
     != solidTractionFvPatchVectorField::typeName
    //!= tractionDisplacementIncrementFvPatchVectorField::typeName
    )
    {
        FatalErrorIn(args.executable())
            << "Bounary condition on " << DU.name()
                <<  " is "
                << DU.boundaryField()[solidPatchID].type()
                << "for fluid -solid interface patch, instead "
            << solidTractionFvPatchVectorField::typeName
            //<< tractionDisplacementIncrementFvPatchVectorField::typeName
                << abort(FatalError);
    }

    //tractionDisplacementIncrementFvPatchVectorField& tForce =
    //      refCast<tractionDisplacementIncrementFvPatchVectorField>
    solidTractionFvPatchVectorField& tForce =
        refCast<solidTractionFvPatchVectorField>
        (
            DU.boundaryField()[solidPatchID]
        );


    // Accumulated fluid interface displacement
    IOobject accumulatedFluidInterfaceDisplacementHeader
    (
        "accumulatedFluidInterfaceDisplacement",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );

    vectorField* accumulatedFluidInterfaceDisplacementPtr = NULL;

    if(accumulatedFluidInterfaceDisplacementHeader.headerOk())
    {
        Pout << "Reading accumulated fluid interface displacement" << endl;

        accumulatedFluidInterfaceDisplacementPtr =
            new vectorIOField
            (
                IOobject
                (
                    "accumulatedFluidInterfaceDisplacement",
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                )
            );
    }
    else
    {
        accumulatedFluidInterfaceDisplacementPtr =
            new vectorIOField
            (
                IOobject
                (
                    "accumulatedFluidInterfaceDisplacement",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                vectorField
                (
                    mesh.boundaryMesh()[fluidPatchID].nPoints(),
                    vector::zero
                )
            );
    }

    vectorField& accumulatedFluidInterfaceDisplacement
        = *accumulatedFluidInterfaceDisplacementPtr;

