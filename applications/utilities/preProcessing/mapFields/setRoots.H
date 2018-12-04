    argList::noParallel();
    argList::validArgs.clear();
    argList::validArgs.append("source dir");

    argList::validOptions.insert("sourceRegion", "name");
    argList::validOptions.insert("targetRegion", "name");

    argList::validOptions.insert("sourceTime", "scalar");
    argList::validOptions.insert("parallelSource", "");
    argList::validOptions.insert("parallelTarget", "");
    argList::validOptions.insert("consistent", "");

    argList args(argc, argv);

    if (!args.check())
    {
         FatalError.exit();
    }

    fileName rootDirTarget(args.rootPath());
    fileName caseDirTarget(args.globalCaseName());

    fileName casePath(args.additionalArgs()[0]);
    fileName rootDirSource = casePath.path();
    fileName caseDirSource = casePath.name();

    word sourceRegionName;
    word targetRegionName;

    args.optionReadIfPresent
    (
        "sourceRegion",
        sourceRegionName,
        fvMesh::defaultRegion
    );

    args.optionReadIfPresent
    (
        "targetRegion",
        targetRegionName,
        fvMesh::defaultRegion
    );

    Info<< "Source: " << rootDirSource << " " << caseDirSource
        << " " << sourceRegionName << nl
        << "Target: " << rootDirTarget << " " << caseDirTarget
        << " " << targetRegionName << endl;

    bool parallelSource = args.optionFound("parallelSource");
    bool parallelTarget = args.optionFound("parallelTarget");
    bool consistent     = args.optionFound("consistent");

