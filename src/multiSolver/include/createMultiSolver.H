    word solverDomain;

    multiSolver multiRun
    (
        multiSolver::multiControlDictName,
        args.rootPath(),
        args.caseName()
    );

    while (multiRun.run())
    { // While brace
        { // dummy brace
