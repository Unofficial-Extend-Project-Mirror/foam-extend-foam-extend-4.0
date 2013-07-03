    Info<< "Create mesh\n" << endl;
    autoPtr<engineTopoChangerMesh> meshPtr
    (
        engineTopoChangerMesh::New
        (
            IOobject
            (
                engineTopoChangerMesh::defaultRegion,
                runTime.timeName(),
                runTime,
                IOobject::MUST_READ
            )
        )
    );

    engineTopoChangerMesh& mesh = meshPtr();

    Switch checkEngineMesh =
    (
        runTime.engineDict().lookup("checkMesh")
    );
