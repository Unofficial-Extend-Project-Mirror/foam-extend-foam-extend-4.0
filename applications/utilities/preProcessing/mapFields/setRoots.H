    argList::noParallel();
    argList::validArgs.clear();
    argList::validArgs.append("source dir");

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

    Info<< "Source: " << rootDirSource << " " << caseDirSource << nl
        << "Target: " << rootDirTarget << " " << caseDirTarget << endl;

    bool parallelSource = false;
    if (args.options().found("parallelSource"))
    {
        parallelSource = true;
    }

    bool parallelTarget = false;
    if (args.options().found("parallelTarget"))
    {
        parallelTarget = true;
    }

    bool consistent = false;
    if (args.options().found("consistent"))
    {
        consistent = true;
    }
