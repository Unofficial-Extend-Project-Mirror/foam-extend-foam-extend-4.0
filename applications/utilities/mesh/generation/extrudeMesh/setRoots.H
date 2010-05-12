    argList::validArgs.clear();
    argList::noParallel();

    argList::validOptions.insert("sourceRoot", "source root");
    argList::validOptions.insert("sourceCase", "source case");
    argList::validOptions.insert("sourcePatch", "source patch");

    argList::validOptions.insert("surface", "surface file");

    argList::validOptions.insert("mergeFaces", "");

    argList args(argc, argv);

    if (!args.check())
    {
         FatalError.exit();
    }

