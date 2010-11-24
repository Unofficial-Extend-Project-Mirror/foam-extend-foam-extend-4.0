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

    if (spacePatchID < 0)
    {
        FatalErrorIn(args.executable())
            << "Space patch not defined. Please make sure that "
            << "the space patch is named as space."
            << abort(FatalError);
    }


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
        if (interface.twoFluids())
        {
            Cb = gSum((1.0 - fluidIndicator.internalField())
               *mesh.C().internalField()*mesh.V())/
                gSum((1.0 - fluidIndicator.internalField())*mesh.V());
        }
        else
        {
            scalar Va = gSum(mesh.V().field());
            vector Ra = 
                gSum(mesh.C().internalField()*mesh.V().field())
               /gSum(mesh.V().field());

            Info << "Ra = " << Ra << endl;
            
            scalar Vb =
              - gSum
                (
                    mesh.Cf().boundaryField()[interface.aPatchID()]
                  & mesh.Sf().boundaryField()[interface.aPatchID()]
                )/3;

            scalar V = Va + Vb;

            Cb = (V*Ra - Va*Ra)/Vb;

            Info << "Current bubble centre: " << Cb << endl;
            
            if (mesh.nGeometricD() != 3)
            {
                FatalErrorIn(args.executable())
                    << "One-fluid bubble centroid calc "
                        << "is not implemented for 2d mesh"
                        << abort(FatalError);
            }
        }

        movingReferenceFrame.add("Cb", Cb);
    }

    // - Bubble centre at the start of calculation
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

        movingReferenceFrame.add("Cbf", Cbf);
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

        movingReferenceFrame.add("XF", XF);
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
    else
    {
        movingReferenceFrame.add("UF", UF);
    }


    // - Bubble acceleration
    dimensionedVector aF ("aF", dimAcceleration, vector::zero);
    if (movingReferenceFrame.found("aF"))
    {
        aF = dimensionedVector
        (
            movingReferenceFrame.lookup("aF")
        );
    }
    else
    {
        movingReferenceFrame.add("aF", aF);
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
        Vb = gSum((1.0 - fluidIndicator.internalField())*mesh.V());

        if (!interface.twoFluids())
        {
            scalar Va = gSum(mesh.V().field());
            
            boundBox box(mesh.C().boundaryField()[spacePatchID]);
            scalar r = (box.max().x() - box.min().x())/2;

            scalar V = 4*M_PI*pow(r,3)/3;

            Vb = V - Va;

            if (mesh.nGeometricD() != 3)
            {
                FatalErrorIn(args.executable())
                    << "One-fluid bubble centroid calc "
                        << "is not implemented for 2d mesh"
                        << abort(FatalError);
            }
        }

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
