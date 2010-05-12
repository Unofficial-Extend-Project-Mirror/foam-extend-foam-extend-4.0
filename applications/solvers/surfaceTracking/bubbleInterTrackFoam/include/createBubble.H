    // Detect the space patch
    label spacePatchID = -1;

    forAll (mesh.boundary(), patchI)
    {
        if (mesh.boundary()[patchI].name() == "space")
        {
            spacePatchID = patchI;

            Info<< "Found space patch. ID: "
                << spacePatchID << endl;
        }
    }


//     if (spacePatchID < 0)
//     {
//         FatalErrorIn(args.executable())
//             << "Space patch not defined. Please make sure that "
//             << "the space patch is named as space."
//             << abort(FatalError);
//     }


    IOdictionary mrfProperties
    (
        IOobject
        (
            "mrfProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    scalar lambdaFf
    (
        readScalar(mrfProperties.lookup("lambdaFf"))
    );

    scalar lambdaF0
    (
        readScalar(mrfProperties.lookup("lambdaF0"))
    );


    // Parameters of moving reference frame
    IOdictionary movingReferenceFrame
    (
        IOobject
        (
            "movingReferenceFrame",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        )
    );

    // - Bubble center
    vector Cb;
    if (movingReferenceFrame.found("Cb"))
    {
        Cb = vector
        (
            movingReferenceFrame.lookup("Cb")
        );
    }
    else
    {
        Cb = gSum((1.0 - interface.fluidIndicator().internalField())*
            mesh.C().internalField()*mesh.V())/
            gSum((1.0 - interface.fluidIndicator().internalField())*mesh.V());
    }

    // - Bubble centre on the start of calculation
    vector Cbf;
    if (movingReferenceFrame.found("Cbf"))
    {
        Cbf = vector
        (
            movingReferenceFrame.lookup("Cbf")
        );
    }
    else
    { 
        Cbf = Cb;
    }

    // - Bubble position vector
    dimensionedVector XF ("XF", dimLength, vector::zero);
    if (movingReferenceFrame.found("XF"))
    {
        XF = dimensionedVector
        (
            movingReferenceFrame.lookup("XF")
        );
    }
    else
    {
        XF = dimensionedVector
        (
            "XF",
            dimLength,
            Cbf
        );
    }


    // - Bubble velocity
    dimensionedVector UF ("UF", dimVelocity, vector::zero);
    if (movingReferenceFrame.found("UF"))
    {
        UF = dimensionedVector
        (
            movingReferenceFrame.lookup("UF")
        );
    }


    // - Bubble acceleration
    dimensionedVector aF ("aF", dimensionSet(0,1,-2,0,0,0,0), vector::zero);
    if (movingReferenceFrame.found("aF"))
    {
        aF = dimensionedVector
        (
            movingReferenceFrame.lookup("aF")
        );
    }


    IOdictionary bubbleProperties
    (
        IOobject
        (
            "bubbleProperties",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        )
    );

    scalar Vb;
    if (bubbleProperties.found("Vb"))
    {
        Vb = readScalar
        (
            bubbleProperties.lookup("Vb")
        );
    }
    else
    {
        Vb = gSum((1.0 - interface.fluidIndicator().internalField())*mesh.V());
        bubbleProperties.add("Vb", Vb);
    }

    scalar Vbf;
    if (bubbleProperties.found("Vbf"))
    {
        Vbf = readScalar
        (
            bubbleProperties.lookup("Vbf")
        );
    }
    else
    { 
        Vbf = Vb;
        bubbleProperties.add("Vbf", Vbf);
    }
