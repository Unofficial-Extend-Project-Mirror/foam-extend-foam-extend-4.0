    OFstream* pathFilePtr = NULL;

    if(Pstream::master())
    { 
        mkDir
        ( 
            args.rootPath()
            /args.globalCaseName()
            /"path"
            /runTime.timeName()
        );
        
        pathFilePtr = 
            new OFstream
            (
                args.rootPath()
                /args.globalCaseName()
                /"path"
                /runTime.timeName()
                /"path"
            );
        
        *pathFilePtr << "Time" << tab;
        
        *pathFilePtr << "Xb" << tab << "Yb" << tab << "Zb" << tab 
            << "Ub" << tab << "ab" << endl;
    }

    OFstream& pathFile = *pathFilePtr;
