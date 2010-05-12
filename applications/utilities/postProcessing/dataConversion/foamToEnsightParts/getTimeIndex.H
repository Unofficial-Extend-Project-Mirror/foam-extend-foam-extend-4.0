// Read time index from */uniform/time,
// but treat 0 and constant specially

    label timeIndex = 0;

    if
    (
        runTime.timeName() != "constant"
     && runTime.timeName() != "0"
    )
    {
        IOobject io
        (
            "time",
            runTime.timeName(),
            "uniform",
            runTime,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        );

        if (io.headerOk())
        {
            IOdictionary timeObject
            (
                IOobject
                (
                    "time",
                    runTime.timeName(),
                    "uniform",
                    runTime,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false
                )
            );

            timeObject.lookup("index") >> timeIndex;
        }
        else
        {
            Info<< "skip ... missing entry " << io.objectPath() << endl;
            continue;
        }
    }

    timeIndices.insert(timeIndex, timeDirs[timeI].value());
    Info<< "\nTime [" << timeIndex << "] = " << runTime.timeName() << nl;

