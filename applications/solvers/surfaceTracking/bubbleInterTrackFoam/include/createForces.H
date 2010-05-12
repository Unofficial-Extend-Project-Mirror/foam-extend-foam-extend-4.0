    OFstream* forcesFilePtr = NULL;

    if(Pstream::master())
    {
        mkDir
        ( 
            args.rootPath()
            /args.globalCaseName()
            /"forces"
            /runTime.timeName()
        );
        
        forcesFilePtr = 
            new OFstream
            (
                args.rootPath()
                /args.globalCaseName()
                /"forces"
                /runTime.timeName()
                /"dragAndLift"
            );

        *forcesFilePtr << "Time" << tab << "drag" << tab << "lift" << endl;
    }

    OFstream& forcesFile = *forcesFilePtr;
