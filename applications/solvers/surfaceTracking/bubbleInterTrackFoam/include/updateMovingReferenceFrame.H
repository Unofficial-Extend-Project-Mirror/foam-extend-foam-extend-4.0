    vector Cb0 = Cb;

    Cb = gSum((1.0 - interface.fluidIndicator().internalField())*
        mesh.C().internalField()*mesh.V())/
        gSum((1.0 - interface.fluidIndicator().internalField())*mesh.V());

    Info << "Current bubble centre: " << Cb << endl;

    dimensionedVector dUF
    (
        "dUF",
        dimVelocity,
        lambdaFf*(Cb-Cbf)/runTime.deltaT0().value()
      + lambdaF0*(Cb-Cb0)/runTime.deltaT0().value()
    );

    XF += (UF + 0.5*dUF)*runTime.deltaT0();
    Info << "Current bubble position: " << XF << endl;


//     dimensionedVector UF00 ("UF00", dimVelocity, UF0.value());
//     UF0 = UF;
    UF += dUF;
    Info << "Current bubble velocity: " << UF << endl;

    aF = dUF/runTime.deltaT0();
//     aF = (3.0*UF - 4.0*UF0 + UF00)/(2.0*runTime.deltaT0());
    Info << "Current bubble acceleration: " << aF << endl;

    if(spacePatchID != -1)
    {
        if
        (
            U.boundaryField()[spacePatchID].type()
         == inletOutletFvPatchVectorField::typeName
        )
        {
            inletOutletFvPatchVectorField& spaceU =
                refCast<inletOutletFvPatchVectorField>
                (
                    U.boundaryField()[spacePatchID]
                );

            spaceU.refValue() = -UF.value();
        }
    }

    // Update movingReferenceFrame dictionary
    if (movingReferenceFrame.found("Cb"))
    {
        movingReferenceFrame.remove("Cb");
    }
    movingReferenceFrame.add("Cb", Cb);

    if (movingReferenceFrame.found("Cbf"))
    {
        movingReferenceFrame.remove("Cbf");
    }
    movingReferenceFrame.add("Cbf", Cbf);

    if (movingReferenceFrame.found("XF"))
    {
        movingReferenceFrame.remove("XF");
    }
    movingReferenceFrame.add("XF", XF);

    if (movingReferenceFrame.found("UF"))
    {
        movingReferenceFrame.remove("UF");
    }
    movingReferenceFrame.add("UF", UF);

    if (movingReferenceFrame.found("aF"))
    {
        movingReferenceFrame.remove("aF");
    }
    movingReferenceFrame.add("aF", aF);
