    Info<< nl << "Create Times" << endl;

    Time runTimeMaster
    (
        Time::controlDictName,
        rootDirMaster,
        caseDirMaster,
        "system",
        "constant",
        !args.optionFound("noFunctionObjects")
    );

    Time runTimeToAdd
    (
        Time::controlDictName,
        rootDirToAdd,
        caseDirToAdd,
        "system",
        "constant",
        !args.optionFound("noFunctionObjects")
    );
