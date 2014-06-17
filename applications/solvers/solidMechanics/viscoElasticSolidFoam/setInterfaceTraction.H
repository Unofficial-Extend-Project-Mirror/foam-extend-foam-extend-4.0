{
    volVectorField totTraction
    (
        IOobject
        (
            "totTraction",
            fluidRunTime.timeName(),
            fluidMesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        fluidMesh
    );


    if
    (
        DU.boundaryField()[interfacePatchID].type()
     != tractionDisplacementIncrementFvPatchVectorField::typeName
    )
    {
        FatalErrorIn(args.executable())
            << "Bounary condition on " << DU.name()
                <<  " is "
                << DU.boundaryField()[interfacePatchID].type()
                << "for patch: "
                << mesh.boundaryMesh()[interfacePatchID].name()
                <<  ", instead "
                << tractionDisplacementIncrementFvPatchVectorField::typeName
                << abort(FatalError);
    }

    tractionDisplacementIncrementFvPatchVectorField& interfacePatchDU =
        refCast<tractionDisplacementIncrementFvPatchVectorField>
        (
            DU.boundaryField()[interfacePatchID]
        );

    interfacePatchDU.traction() =
        interpolator.faceInterpolate
        (
            totTraction.boundaryField()[fluidInterfacePatchID]
        );
}


