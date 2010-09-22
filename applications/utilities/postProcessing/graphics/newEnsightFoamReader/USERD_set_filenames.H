//======================================================================
//  Setting filenames
//======================================================================
int USERD_set_filenames
(
    char filename_1[],
    char filename_2[],
    char the_path[],
    int swapbytes
)
{
#ifdef ENSIGHTDEBUG
    Info << "Entering: USERD_set_filenames" << endl << flush;
#endif

    char tmp[100];

    label lRoot = strlen(the_path);
    label lCase = strlen(filename_1);

    bool cleared = false;

    while (!cleared)
    {
        lRoot = strlen(the_path);
        lCase = strlen(filename_1);

        // remove the last '/' from rootDir
        if (the_path[lRoot-1] == '/')
        {
            the_path[lRoot-1] = (char)NULL;
        }
        else
        {
            cleared = true;
        }
    }

    rootDir = the_path;

    // the path is pre-pended to filename_1
    // 1 is the 'Geometry' : 2 the 'Result' which is null here
    // since two_field is FALSE
    for (label i=0; i<lCase-lRoot;i++)
    {
        tmp[i] = filename_1[i+1+lRoot];
    }
    caseDir = tmp;

    if (!isDir(rootDir/caseDir))
    {
       Info<< rootDir/caseDir << " is not a valid directory."
           << endl;
       return Z_ERR;
    }


    // construct the global pointers to the database and mesh

    delete meshPtr;
    delete runTimePtr;

    runTimePtr = new Time
    (
        Time::controlDictName,
        rootDir,
        caseDir
    );

    Time& runTime = *runTimePtr;

    meshPtr = new fvMesh
    (
        IOobject
        (
            fvMesh::defaultRegion,
            runTime.timeName(),
            runTime
        )
    );

    // set the available number of time-steps
    TimeList = (const instantList&)Foam::Time::findTimes(rootDir/caseDir);

    Num_time_steps = TimeList.size() - 1;

    nPatches = meshPtr->boundaryMesh().size();

    // set the number of fields and store their names
    // a valid field must exist for all time-steps
    runTime.setTime(TimeList[TimeList.size()-1], TimeList.size()-1);
    IOobjectList objects(*meshPtr, runTime.timeName());

    fieldNames = (const wordList&)objects.names();

    // because of the spray being a 'field' ...
    // get the availabe number of variables and
    // check for type (scalar/vector/tensor)

    label nVar = 0;
    wordList scalars = objects.names(scalarName);

    for (label n=0; n<fieldNames.size(); n++)
    {
        bool isitScalar = false;
        forAll(scalars,i)
        {
            if (fieldNames[n] == scalars[i])
            {
                isitScalar = true;
                var2field[nVar++] = n;
            }
        }
        isScalar[n] = isitScalar;
    }

    wordList vectors = objects.names(vectorName);

    for (label n=0; n<fieldNames.size(); n++)
    {
        bool isitVector = false;
        forAll(vectors,i)
        {
            if (fieldNames[n] == vectors[i])
            {
                isitVector = true;
                var2field[nVar++] = n;
            }
        }
        isVector[n] = isitVector;
    }

    wordList tensors = objects.names(tensorName);

    for (label n=0; n<fieldNames.size(); n++)
    {
        bool isitTensor = false;
        forAll(tensors,i)
        {
            if (fieldNames[n] == tensors[i])
            {
                isitTensor = true;
                var2field[nVar++] = n;
            }
        }
        isTensor[n] = isitTensor;
    }
    
    bool lagrangianNamesFound = false;
    label n = 0;
    while ((!lagrangianNamesFound) && (n<Num_time_steps))
    {
        runTime.setTime(TimeList[n+1], n+1);

        Cloud<passiveParticle> lagrangian(*meshPtr);

        n++;
        if (lagrangian.size()>0)
        {
            lagrangianNamesFound = true;
        }
    }

    IOobject sprayHeader
    (
        "positions",
        runTime.timeName(),
        "lagrangian",
        runTime,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false
    );


    if (sprayHeader.headerOk())
    {
        Info << "[Found lagrangian]" << endl;

        delete sprayPtr;

        sprayPtr = new Cloud<passiveParticle>(*meshPtr);

        IOobjectList objects(*meshPtr, runTime.timeName(), "lagrangian");

        lagrangianScalarNames =
            (const wordList&)objects.names(sprayScalarFieldName);
        lagrangianVectorNames =
            (const wordList&)objects.names(sprayVectorFieldName);

        isSpray[fieldNames.size()] = true;

        nSprayVariables += lagrangianScalarNames.size();
        nSprayVariables += lagrangianVectorNames.size();

        Num_unstructured_parts++;
    }

    Current_time_step = Num_time_steps;
    runTime.setTime(TimeList[Current_time_step], Current_time_step);

    Num_variables = nVar + nSprayVariables;


    // Check if second mesh (region) exist. 
    IOobject ensightReaderDictHeader
    (
        IOobject
        (
            "ensightReaderDict",
            runTime.constant(),
            *meshPtr,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    word secondRegionName("solid");

    if(ensightReaderDictHeader.headerOk())
    {
        IOdictionary ensightReaderDict
        (
            IOobject
            (
                "ensightReaderDict",
                runTime.constant(),
                *meshPtr,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        if (ensightReaderDict.found("secondRegionName"))
        {
            secondRegionName = word
            (
                ensightReaderDict.lookup("secondRegionName")
            );
        }

        if (ensightReaderDict.found("scalingFactor"))
        {
            scalingFactor = readScalar
            (
                ensightReaderDict.lookup("scalingFactor")
            );
        }
    }

    objectRegistry secondMeshObjReg
    (
        IOobject
        (
            secondRegionName,
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ
        )
    );

    IOobject pointsHeader
    (
        "points",
        runTime.constant(),
        polyMesh::meshSubDir,
        secondMeshObjReg
    );

    if(pointsHeader.headerOk())
    {
        secondMeshPtr = new fvMesh
        (
            IOobject
            (
                secondRegionName,
                runTime.timeName(),
                runTime,
                IOobject::MUST_READ
            )
        );        
    }

    if(secondMeshPtr)
    {
        Num_unstructured_parts += 1;

        nSecondMeshPatches = secondMeshPtr->boundaryMesh().size();

        IOobjectList secondMeshObjects(*secondMeshPtr, runTime.timeName());

        secondMeshFieldNames = (const wordList&)secondMeshObjects.names();

        label nVar = 0;
        wordList scalars = secondMeshObjects.names(scalarName);

        for (label n=0; n<secondMeshFieldNames.size(); n++)
        {
            bool isitScalar = false;
            forAll(scalars,i)
            {
                if (secondMeshFieldNames[n] == scalars[i])
                {
                    isitScalar = true;
                    secondMeshVar2field[nVar++] = n;
                }
            }
            secondMeshIsScalar[n] = isitScalar;
        }

        wordList vectors = secondMeshObjects.names(vectorName);

        for (label n=0; n<secondMeshFieldNames.size(); n++)
        {
            bool isitVector = false;
            forAll(vectors,i)
            {
                if (secondMeshFieldNames[n] == vectors[i])
                {
                    isitVector = true;
                    secondMeshVar2field[nVar++] = n;
                }
            }
            secondMeshIsVector[n] = isitVector;
        }

        wordList tensors = secondMeshObjects.names(tensorName);

        for (label n=0; n<secondMeshFieldNames.size(); n++)
        {
            bool isitTensor = false;
            forAll(tensors,i)
            {
                if (secondMeshFieldNames[n] == tensors[i])
                {
                    isitTensor = true;
                    secondMeshVar2field[nVar++] = n;
                }
            }
            secondMeshIsTensor[n] = isitTensor;
        }

        Num_variables += nVar;
        nSecondMeshVariables = nVar;
    }

    secondMeshPartNum = nPatches+2;

    if(sprayPtr)
    {
        secondMeshPartNum += 1;
    }


    // Check if finite area mesh exists
    IOobject faMeshBoundaryIOobj
    (
        "boundary",
        meshPtr->time().findInstance
        (
            meshPtr->dbDir()/fvMesh::meshSubDir, "boundary"
        ),
        faMesh::meshSubDir,
        *meshPtr,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );


    if(faMeshBoundaryIOobj.headerOk())
    {
        Info << "\nFound finite area mesh" << endl;
        faMeshPtr = new faMesh(*meshPtr);
    }


    if(faMeshPtr)
    {
        Num_unstructured_parts += 1;

        IOobjectList faMeshObjects(faMeshPtr->time(), runTime.timeName());

        faMeshFieldNames = (const wordList&)faMeshObjects.names();

        Info << faMeshFieldNames << endl;

        label nVar = 0;
        wordList scalars = faMeshObjects.names(faScalarName);

        for (label n=0; n<faMeshFieldNames.size(); n++)
        {
            bool isitScalar = false;
            forAll(scalars,i)
            {
                if (faMeshFieldNames[n] == scalars[i])
                {
                    isitScalar = true;
                    faMeshVar2field[nVar++] = n;
                }
            }
            faMeshIsScalar[n] = isitScalar;
        }

        wordList vectors = faMeshObjects.names(faVectorName);

        for (label n=0; n<faMeshFieldNames.size(); n++)
        {
            bool isitVector = false;
            forAll(vectors,i)
            {
                if (faMeshFieldNames[n] == vectors[i])
                {
                    isitVector = true;
                    faMeshVar2field[nVar++] = n;
                }
            }
            faMeshIsVector[n] = isitVector;
        }

        wordList tensors = faMeshObjects.names(faTensorName);

        for (label n=0; n<faMeshFieldNames.size(); n++)
        {
            bool isitTensor = false;
            forAll(tensors,i)
            {
                if (faMeshFieldNames[n] == tensors[i])
                {
                    isitTensor = true;
                    faMeshVar2field[nVar++] = n;
                }
            }
            faMeshIsTensor[n] = isitTensor;
        }

        Num_variables += nVar;
        nFaMeshVariables = nVar;
    }


    Numparts_available = Num_unstructured_parts + Num_structured_parts 
        + nPatches + nSecondMeshPatches;


#ifdef ENSIGHTDEBUG
    Info << "Leaving: USERD_set_filenames" << endl << flush;
#endif

    return Z_OK;
}

