
    fileName casePath(args.additionalArgs()[0]);
    fileName fluidRootDir = casePath.path();
    fileName fluidCaseDir = casePath.name();

    Info<< "Fluid: " << fluidRootDir << " " << fluidCaseDir << endl;

    Time fluidRunTime
    (
        Time::controlDictName,
        fluidRootDir,
        fluidCaseDir
    );

    fluidRunTime.setDeltaT(runTime.deltaT().value());

    fvMesh fluidMesh
    (
        IOobject
        (
            fvMesh::defaultRegion,
            fluidRunTime.timeName(),
            fluidRunTime
        )
    );

    label fluidInterfacePatchID =
        fluidMesh.boundaryMesh().findPatchID("interface");

    if(fluidInterfacePatchID == -1)
    {
        FatalErrorIn(args.executable())
            << "Can't find patch with name interface" << abort(FatalError);
    }

    patchToPatchInterpolation interpolator
    (
        fluidMesh.boundaryMesh()[fluidInterfacePatchID],
        mesh.boundaryMesh()[interfacePatchID],
        intersection::VISIBLE
    );

    scalarField error =
        mag
        (
            mesh.C().boundaryField()[interfacePatchID]
          - interpolator.faceInterpolate(fluidMesh.C().boundaryField()[fluidInterfacePatchID])
        );

    Info << "Patch-to-patch interpolation error, max: " << max(error) << ", avg: " << average(error) << endl;
