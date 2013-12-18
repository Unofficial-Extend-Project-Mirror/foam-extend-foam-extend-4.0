
    argList::noParallel();
    argList::validArgs.clear();
//     argList::validArgs.append("fluid case dir");

    argList args(argc, argv);

    if (!args.checkRootCase())
    {
         FatalError.exit();
    }
