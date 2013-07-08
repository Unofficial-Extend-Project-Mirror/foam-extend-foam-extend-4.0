/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "multiSolver.H"
#include "tuple2Lists.H"
#include "OFstream.H"
#include "Pstream.H"

// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(multiSolver, 0);
}


template<>
const char* Foam::NamedEnum<Foam::multiSolver::initialStartFromControls, 9>::names[] =
{
    "firstTime",
    "firstTimeInStartDomain",
    "firstTimeInStartDomainInStartSuperLoop",
    "startTime",
    "startTimeInStartDomain",
    "startTimeInStartDomainInStartSuperLoop",
    "latestTime",
    "latestTimeInStartDomain",
    "latestTimeInStartDomainInStartSuperLoop"
};

const Foam::NamedEnum<Foam::multiSolver::initialStartFromControls, 9>
    Foam::multiSolver::initialStartFromControlsNames_;


template<>
const char* Foam::NamedEnum<Foam::multiSolver::finalStopAtControls, 7>::names[] =
{
    "endTime",
    "endTimeInEndDomain",
    "endTimeInEndDomainInEndSuperLoop",
    "endSuperLoop",
    "writeNow",
    "noWriteNow",
    "nextWrite"
};

const Foam::NamedEnum<Foam::multiSolver::finalStopAtControls, 7>
    Foam::multiSolver::finalStopAtControlsNames_;


template<>
const char* Foam::NamedEnum<Foam::multiSolver::startFromControls, 4>::names[] =
{
    "firstTime",
    "startTime",
    "latestTimeThisDomain",
    "latestTimeAllDomains"
};

const Foam::NamedEnum<Foam::multiSolver::startFromControls, 4>
    Foam::multiSolver::startFromControlsNames_;


template<>
const char* Foam::NamedEnum<Foam::multiSolver::stopAtControls, 7>::names[] =
{
    "endTime",
    "noWriteNow",
    "writeNow",
    "nextWrite",
    "iterations",
    "solverSignal",
    "elapsedTime"
};

const Foam::NamedEnum<Foam::multiSolver::stopAtControls, 7>
    Foam::multiSolver::stopAtControlsNames_;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::word Foam::multiSolver::multiControlDictName("multiControlDict");


void Foam::multiSolver::setUpParallel()
{
    if (Pstream::master())
    {
        fileNameList roots(Pstream::nProcs());
        
        roots[0] = multiDictRegistry_.rootPath();
        manageLocalRoot_ = true;
        
        // Receive from slaves
        for
        (
            int slave=Pstream::firstSlave();
            slave<=Pstream::lastSlave();
            slave++
        )
        {
            IPstream fromSlave(Pstream::blocking, slave);
            roots[slave] = fileName(fromSlave);
        }
        
        // Distribute
        for
        (
            int slave=Pstream::firstSlave();
            slave<=Pstream::lastSlave();
            slave++
        )
        {
            OPstream toSlave(Pstream::blocking, slave);
            if (roots[slave] != roots[slave - 1])
            {
                toSlave << true;
            }
            else
            {
                toSlave << false;
            }
        }
    }
    else
    {
        // Send to master
        {
            OPstream toMaster(Pstream::blocking, Pstream::masterNo());
            toMaster << fileName(multiDictRegistry_.rootPath());
        }
        // Receive from master
        {
            IPstream fromMaster(Pstream::blocking, Pstream::masterNo());
            manageLocalRoot_ = readBool(fromMaster);
        }
    }
}


void Foam::multiSolver::synchronizeParallel() const
{
    if (Pstream::master())
    {
        // Give go signal
        for
        (
            int slave=Pstream::firstSlave();
            slave<=Pstream::lastSlave();
            slave++
        )
        {
            OPstream toSlave(Pstream::blocking, slave);
            toSlave << true;
        }
    }
    else
    {
        // Recieve go signal
        {
            IPstream fromMaster(Pstream::blocking, Pstream::masterNo());
            bool okayToGo(readBool(fromMaster));
        }
    }
}


Foam::word Foam::multiSolver::setLocalEndTime()
{
    word stopAtSetting("endTime");
    switch (stopAt_)
    {
        case msaEndTime:
            // do nothing
            break;
        case msaNoWriteNow:
            stopAtSetting = "noWriteNow";
            finished_ = true;
            break;
        case msaWriteNow:
            stopAtSetting = "writeNow";
            finished_ = true;
            break;
        case msaNextWrite:
            stopAtSetting = "nextWrite";
            finished_ = true;
            break;
        case msaIterations:
            endTime_ = deltaT_ * iterations_ + startTime_;
            break;
        case msaSolverSignal:
            endTime_ = VGREAT;
            break;
        case msaElapsedTime:
            endTime_ = startTime_ + elapsedTime_;
            break;
    }

    // Modify endTime_ if it exceeds finalEndTime
    switch (finalStopAt_)
    {
        case mfsEndTime:
            if ((endTime_ + globalTimeOffset_) >= finalEndTime_)
            {
                endTime_ = finalEndTime_ - globalTimeOffset_;
                finished_ = true;
                if ((startTime_ + globalTimeOffset_) >= finalEndTime_)
                {
                    // Initialized beyond end
                    stopAtSetting = "noWriteNow";
                }
            }

            break;
        case mfsEndTimeInEndDomain:
            if 
            (
                (currentSolverDomain_ == endDomain_)
             && (endTime_ >= finalEndTime_)
            )
            {
                endTime_ = finalEndTime_;
                finished_ = true;
                if (startTime_ >= finalEndTime_)
                {
                    // Initialized beyond end
                    stopAtSetting = "noWriteNow";
                }
            }
            break;
        case mfsEndTimeInEndDomainInEndSuperLoop:
            if (currentSolverDomain_ == endDomain_)
            {
                if
                (
                    (superLoop_ >= endSuperLoop_)
                 && (endTime_ >= finalEndTime_)
                )
                {
                    endTime_ = finalEndTime_;
                    finished_ = true;
                    if (startTime_ > finalEndTime_)
                    {
                        // Initialized beyond end
                        stopAtSetting = "noWriteNow";
                    }
                }
            }
            // Cannot check for beyond end initialization with superLoops
            // because operator++ allows the superLoop to increment by more
            // than 1
            break;
        case mfsEndSuperLoop:
            if (superLoop_ > endSuperLoop_)
            {
                stopAtSetting = "noWriteNow";
                finished_ = true;
            }
            break;
        case mfsWriteNow:
            finished_ = true;
            stopAtSetting = "writeNow";
            break;
        case mfsNoWriteNow:
            finished_ = true;
            stopAtSetting = "noWriteNow";
            break;
        case mfsNextWrite:
            stopAtSetting = "nextWrite";
            finished_ = true;
            break;
    }
    return stopAtSetting;
}
    
void Foam::multiSolver::checkTimeDirectories() const
{
    forAll(prefixes_, i)
    {
        if (prefixes_[i] == "default") continue;
        if ((prefixes_[i] == "all") || (prefixes_[i] == "root"))
        {
            FatalErrorIn("multiSolver::checkTimeDirectories")
                << "'all' or 'root' solverDomain name found in "
                << "multiControlDict.  These two names are prohibitted."
                << abort(FatalError);
        }
        fileName checkMe
        (
            multiDictRegistry_.path()/"multiSolver"/prefixes_[i]/"initial/0"
        );
        if
        (
            !exists(checkMe)
        )
        {
            FatalErrorIn("multiSolver::checkTimeDirectories")
                << "Initial time directory missing for solver domain ["
                << prefixes_[i] << "].  Expecting " << checkMe
                << abort(FatalError);
        }
    }
}


void Foam::multiSolver::swapDictionaries(const word& solverDomainName)
{
    forAll(multiDicts_, i)
    {
        IOdictionary newMultiDict
        (
            IOobject
            (
                multiDicts_[i].lookup("dictionaryName"),
                multiDicts_[i].instance(),
                multiDicts_[i].local(),
                multiDictRegistry_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                false
            )
        );
        if (multiDicts_[i].subDict("multiSolver").found("default"))
        {
            newMultiDict.merge
            (
                multiDicts_[i].subDict("multiSolver").subDict("default")
            );
        }
        
        if (multiDicts_[i].subDict("multiSolver").found(solverDomainName))
        {
            if (multiDicts_[i].subDict("multiSolver")
                .subDict(solverDomainName).found("sameAs"))
            {
                word sameAsSolverDomain(multiDicts_[i].subDict("multiSolver")
                    .subDict(solverDomainName).lookup("sameAs"));
                if (multiControlDict_.subDict("solverDomains")
                    .found(sameAsSolverDomain))
                {
                    newMultiDict.merge
                    (
                        multiDicts_[i].subDict("multiSolver")
                            .subDict(sameAsSolverDomain)
                    );
                }
                else
                {
                    FatalIOErrorIn
                    (
                        "multiSolver::swapDictionaries", multiDicts_[i]
                    )
                        << "'sameAs' solverDomain name " << sameAsSolverDomain
                        << " not found."
                        << exit(FatalIOError);
                }
            }
            else
            newMultiDict.merge
            (
                multiDicts_[i].subDict("multiSolver").subDict(solverDomainName)
            );
        }

        newMultiDict.regIOobject::write();
    }
}


void Foam::multiSolver::swapBoundaryConditions
(
    const fileName& dataSourcePath,
    const word& intoSolverDomain
)
{
    fileName bcFilePath
    (
        multiDictRegistry_.path()/"multiSolver"/intoSolverDomain/"initial/0"
    );

    fileName dataSourceInitialConditions
    (
        multiDictRegistry_.path()/"multiSolver"/currentSolverDomain_
            /"initial/0"
    );

    instantList il(Time::findTimes(dataSourcePath.path()));

    fileName firstDataSourcePath
    (
        dataSourcePath.path()/il[Time::findClosestTimeIndex(il,-1.0)].name()
    );

    fileNameList dirEntries
    (
        readDir(dataSourcePath, fileName::FILE)
    );
    
    word headerClassName;
    forAll(dirEntries, i)
    {
        // Ignore this file if it isn't in case/prefix/initial/0
        if (!exists(bcFilePath/dirEntries[i])) continue;

        IFstream bc(bcFilePath/dirEntries[i]);
        IFstream data(dataSourcePath/dirEntries[i]);

        // Find the headerClassName for pseudo regIOobject::write later
        while (!bc.eof())
        {
            token nextToken(bc);
            if (nextToken.isWord() && nextToken.wordToken() == "class")
            {
                break;
            }
        }
        headerClassName = word(token(bc).wordToken());
        bc.rewind();

        dictionary bcDict(bc);
        dictionary dataDict(data);
        
        if
        (
            !bcDict.found("dimensions")
         || !bcDict.found("internalField")
         || !bcDict.found("boundaryField")
         || !dataDict.found("dimensions")
         || !dataDict.found("internalField")
         || !dataDict.found("boundaryField")
        )
        {
            // Data source or BC source not a proper geometricField file
            continue;
        }

        dimensionSet bcDims(bcDict.lookup("dimensions"));
        dimensionSet dataDims(dataDict.lookup("dimensions"));

        if ((bcDims != dataDims) && dimensionSet::debug)
        {
            FatalErrorIn("multiSolver::swapBoundaryConditions")
            << "Dimensions do not match in geometricFields with the same "
            << "name.  Solver domain [" << intoSolverDomain << "] has "
            << bcDims << " and the previous domain has " << dataDims << "."
            << abort(FatalError);
        }
        
        dictionary outputDict(bcDict);

        outputDict.set("internalField", dataDict.lookup("internalField"));

        wordList dataPatches(dataDict.subDict("boundaryField").toc());
        wordList bcPatches(bcDict.subDict("boundaryField").toc());
        sort(dataPatches);
        sort(bcPatches);
        
        if (dataPatches.size() != bcPatches.size())
        {
            FatalErrorIn("multiSolver::swapBoundaryConditions")
            << "Boundary fields do not match.  Solver domain [" 
            << intoSolverDomain << "] has " << bcPatches.size() << " patches "
            << "and the previous domain has " << dataPatches.size() << "."
            << abort(FatalError);
        }

        forAll(dataPatches, j)
        {
            if (dataPatches[j] != bcPatches[j])
            {
                FatalErrorIn("multiSolver::swapBoundaryConditions")
                << "Boundary fields do not match.  Solver domain [" 
                << intoSolverDomain << "] has:" << bcPatches << " patches "
                << "and the previous domain has:" << dataPatches << "."
                << abort(FatalError);
            }
            if (exists(firstDataSourcePath/dirEntries[i]))
            {
                IFstream firstDataStream(firstDataSourcePath/dirEntries[i]);
                dictionary firstDict(firstDataStream);
                
                // Check for 'multiSolverRemembering' entries, copy them from
                // the earliest time (this superLoop) to the outputDict
                if
                (
                    firstDict.subDict("boundaryField")
                        .subDict(bcPatches[j]).found("multiSolverRemembering")
                )
                {
                    wordList msr
                    (
                        firstDict.subDict("boundaryField")
                            .subDict(bcPatches[j])
                            .lookup("multiSolverRemembering")
                    );
                    forAll(msr, k)
                    {
                        if (!firstDict.found(msr[k]))
                        {
                            FatalIOErrorIn
                            (
                                "multiSolver::swapBoundaryConditions",
                                firstDict
                            )
                                << "'multiSolverRemember' word '" << msr[k]
                                << "' missing from boundary patch '"
                                << dataPatches[j] << "' while switching to "
                                << currentSolverDomain_ << ".  This may be "
                                << "the result of manual editting datafiles "
                                << "or data corruption.  If the problem "
                                << "persists, this is a bug."
                                << exit(FatalIOError);
                        }

                        outputDict.subDict("boundaryField")
                            .subDict(bcPatches[j]).set
                        (
                            msr[k],
                            firstDict.subDict("boundaryField")
                                .subDict(bcPatches[j]).lookup(msr[k])
                        );
                    }
                    outputDict.subDict("boundaryField")
                        .subDict(bcPatches[j]).set
                    (
                        "multiSolverRemembering",
                        msr
                    );
                }

                // Check for "multiSolverRemember" fields, copy them from
                // latestTime to outputDict, append their names to multiSolver-
                // Remembering
                if (exists(dataSourceInitialConditions/dirEntries[i]))
                {
                    IFstream ic(dataSourceInitialConditions/dirEntries[i]);
                    dictionary icDict(ic);
                    if
                    (
                        icDict.subDict("boundaryField")
                            .subDict(bcPatches[j]).found("multiSolverRemember")
                    )
                    {
                        wordList remember
                        (
                            icDict
                                .subDict("boundaryField")
                                .subDict(bcPatches[j])
                                .lookup("multiSolverRemember")
                        );

                        forAll(remember, k)
                        {
                            if (!dataDict.found(remember[k]))
                            {
                                FatalIOErrorIn
                                (
                                    "multiSolver::swapBoundaryConditions",
                                    dataDict
                                )
                                    << "'multiSolverRemember' wordList found, "
                                    << "but keyword '" << remember[k] << "' not"
                                    << "present in dictionary for "
                                    << dirEntries[i]
                                    << exit(FatalIOError);
                            }
                            
                            outputDict
                                .subDict("boundaryField")
                                .subDict(bcPatches[j]).set
                            (
                                remember[j],
                                dataDict.subDict("boundaryField")
                                    .subDict(bcPatches[j]).lookup(remember[k])
                            );
                        }

                        wordList remembering(remember);

                        if
                        (
                            firstDict.subDict("boundaryField")
                                .subDict(bcPatches[j])
                                .found("multiSolverRemembering")
                        )
                        {
                            wordList msr
                            (
                                firstDict.subDict("boundaryField")
                                    .subDict(bcPatches[j])
                                    .lookup("multiSolverRemembering")
                            );
                            remembering.setSize(remember.size() + msr.size());
                            forAll(msr, l)
                            {
                                remembering[remember.size() + l] = msr[l];
                            }
                        }
                        
                        outputDict
                            .subDict("boundaryField")
                            .subDict(bcPatches[j]).set
                        (
                            "multiSolverRemembering",
                            remembering
                        );
                    }
                }
            } // End multiSolverRemember implementation
        } // end cycle through patches

        // Here we are cheating a regIOobject::write from a non regIOobject.
        // This allows us to change the header as we want. * High maintenance*
        OFstream os
        (
            multiDictRegistry_.path()/
                multiDictRegistry_.timeName()/
                dirEntries[i]
        );
        IOobject::writeBanner(os);
        os  << "FoamFile\n{\n"
            << "    version     " << os.version() << ";\n"
            << "    format      " << os.format() << ";\n"
            << "    class       " << headerClassName << ";\n";

        os  << "    object      " << dirEntries[i] << ";\n"
            << "}" << nl;

        IOobject::writeDivider(os);
        os  << endl;
        outputDict.write(os);
    } // end cycle through files
}


void Foam::multiSolver::readAllMultiDicts()
{
    fileName systemPath
    (
        multiDictRegistry_.path()/multiDictRegistry_.system()
    );
    fileName constantPath
    (
        multiDictRegistry_.path()/multiDictRegistry_.constant()
    );
    
    label done(0);
    while (done <= 0)
    {
        readMultiDictDirectory(systemPath);
        readMultiDictDirectory(constantPath);

        // Sub directories under system
        fileNameList dirEntries
        (
            readDir(systemPath, fileName::DIRECTORY)
        );

        forAll(dirEntries, i)
        {
            readMultiDictDirectory
            (
                systemPath/dirEntries[i],
                dirEntries[i]
            );
        }

        // Sub directories under constant
        dirEntries = readDir
        (
            constantPath,
            fileName::DIRECTORY
        );
        
        forAll(dirEntries, i)
        {
            readMultiDictDirectory
            (
                constantPath/dirEntries[i],
                dirEntries[i]
            );
        }
    
        // Add local root for parallel runs and repeat
        if (manageLocalRoot_)
        {
            if (done == 0)
            {
                systemPath =
                    multiDictRegistry_.path().path()
                    /multiDictRegistry_.system();
                constantPath =
                    multiDictRegistry_.path().path()
                    /multiDictRegistry_.constant();
                done = -1;
            }
            else
            {
                done = 1;
            }
        }
        else
        {
            done = 1;
        }
    }
}


void Foam::multiSolver::readMultiDictDirectory
(
    const fileName& sourcePath,
    const word& local
)
{
    fileNameList dirEntries(readDir(sourcePath, fileName::FILE));
    forAll(dirEntries, i)
    {
        if
        (
            (dirEntries[i](5) == "multi")
         && (dirEntries[i] != multiControlDictName)
         && (dirEntries[i] != "multiSolverTime")
        )
        {
            IFstream is(sourcePath/dirEntries[i]);
            dictionary candidate(is);
            
            if
            (
                candidate.found("dictionaryName")
             && candidate.found("multiSolver")
            )
            {
                multiDicts_.setSize(multiDicts_.size() + 1);
                multiDicts_.set
                (
                    multiDicts_.size() - 1,
                    new IOdictionary
                    (
                        IOobject
                        (
                            dirEntries[i],
                            sourcePath.name(),
                            local,
                            multiDictRegistry_,
                            IOobject::MUST_READ,
                            IOobject::NO_WRITE
                        )
                    )
                );
            }
        }
    }
}


void Foam::multiSolver::readIfModified()
{
    if (multiDictsRunTimeModifiable_)
    {
        multiDictRegistry_.readModifiedObjects();
        setMultiSolverControls();
    }
}


Foam::timeCluster Foam::multiSolver::parseConditionedFile
(
    const word& pcFile,
    const instant& inst
) const
{
    // solverDomain@superLoop@globalOffset@preConName
#   ifdef FULLDEBUG
    if (!pcFile.size())
    {
        FatalErrorIn("multiSolver::parseConditionedFile")
            << "Empty preConditioned fileName."
            << abort(FatalError);
    }
#   endif

    // Find first @
    string::size_type first = pcFile.find("@");
#   ifdef FULLDEBUG
    if (first == string::npos)
    {
        FatalErrorIn("multiSolver::parseConditionedFile")
            << "Bad preConditioned fileName: " << pcFile
            << abort(FatalError);
    }
#   endif

    // Find second @
    string::size_type second = pcFile(first + 1, pcFile.size() - first - 1)
        .find("@");
#   ifdef FULLDEBUG
    if (second == string::npos)
    {
        FatalErrorIn("multiSolver::parseConditionedFile")
            << "Bad preConditioned fileName: " << pcFile
            << abort(FatalError);
    }
#   endif

    // Find third @
    string::size_type third = pcFile
    (
        first + second + 2,
        pcFile.size() - first - second - 2
    ).find("@");

#   ifdef FULLDEBUG
    if
    (
        third == string::npos
     || pcFile.size() == first + second + third + 3
    )
    {
        FatalErrorIn("multiSolver::parseConditionedFile")
            << "Bad preConditioned fileName: " << pcFile
            << abort(FatalError);
    }
#   endif

    word solverDomain
    (
        pcFile(first)
    );
    
    IStringStream superLoopStream(pcFile(first + 1, second));
    token superLoopToken(superLoopStream);

#   ifdef FULLDEBUG
    if (!superLoopToken.isLabel() || !superLoopStream.eof())
    {
        FatalErrorIn("multiSolver::parseConditionedFile")
            << "Bad preConditioned fileName: " << pcFile
            << abort(FatalError);
    }
#   endif

    label superLoop
    (
        superLoopToken.labelToken()
    );
    
    IStringStream globalOffsetStream(pcFile(first + second + 2, third));
    token globalOffsetToken(globalOffsetStream);

#   ifdef FULLDEBUG
    if (!globalOffsetToken.isNumber() || !globalOffsetStream.eof())
    {
        FatalErrorIn("multiSolver::parseConditionedFile")
            << "Bad preConditioned fileName: " << pcFile
            << abort(FatalError);
    }
#   endif

    scalar globalOffset
    (
        globalOffsetToken.number()
    );
    
    word preConName
    (
        pcFile
        (
            first + second + third + 3,
            pcFile.size() - first - second - third - 3
        )
    );
    
    return timeCluster
    (
        instantList(1, inst),
        globalOffset,
        superLoop,
        solverDomain,
        preConName
    );
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::multiSolver::multiSolver
(
    const dictionary& dict,
    const fileName& rootPath,
    const fileName& caseName,
    const word& systemName,
    const word& constantName
)
:
    dcd_(dict),
    
    multiDictRegistry_
    (
        dcd_,
        rootPath,
        caseName,
        systemName,
        constantName
    ),
    
    multiControlDict_
    (
        IOobject
        (
            multiControlDictName,
            multiDictRegistry_.system(),
            multiDictRegistry_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        ),
        dict
    ),
#include "multiSolverInit.H"
{
    if (Pstream::parRun())
    {
        setUpParallel();
    }
    readAllMultiDicts();
    checkTimeDirectories();
    setMultiSolverControls();
}


Foam::multiSolver::multiSolver
(
    const word& multiControlDictName,
    const fileName& rootPath,
    const fileName& caseName,
    const word& systemName,
    const word& constantName
)
:
    dcd_(rootPath/caseName/systemName/multiControlDictName),
    
    multiDictRegistry_
    (
        dcd_,
        rootPath,
        caseName,
        systemName,
        constantName
    ),

    multiControlDict_
    (
        IOobject
        (
            multiControlDictName,
            multiDictRegistry_.system(),
            multiDictRegistry_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    ),
#include "multiSolverInit.H"
{
    if (Pstream::parRun())
    {
        setUpParallel();
    }
    readAllMultiDicts();
    checkTimeDirectories();
    setMultiSolverControls();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::multiSolver::~multiSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::multiSolver::preCondition(const word& processor)
{
    fileName path(multiDictRegistry_.path());
    if (processor.size())
    {
        path = path/processor;
    }

    // Remove existing time directories from root
    purgeTimeDirs(path);

    //Read all data in path/multiSolver
    timeClusterList tclSource
    (
        readAllTimes(processor)
    );
    tclSource.purgeEmpties();
    forAll(tclSource, tc)
    {
        forAll(tclSource[tc], inst)
        {        
            fileName sourcePath;
            fileName destPath;

            // Destination path goes to -1 for initial conditions because
            // reconstructPar omits 0
            if (tclSource[tc].superLoop() == -1)
            {
                sourcePath = path/"multiSolver"
                    /tclSource[tc].solverDomainName()
                    /"initial/0";
                destPath = path/"-1";
            }
            else
            {
                sourcePath = path/"multiSolver"
                    /tclSource[tc].solverDomainName()
                    /name(tclSource[tc].superLoop())
                    /tclSource[tc][inst].name();
                destPath = path/tclSource[tc][inst].name();
            }
            mkDir(destPath);
            
            fileNameList rootFiles
            (
                readDir(sourcePath, fileName::FILE)
            );
            forAll(rootFiles, rf)
            {
                cp
                (
                    sourcePath/rootFiles[rf],
                    destPath/tclSource[tc].solverDomainName()
                  + "@" + name(tclSource[tc].superLoop())
                  + "@" + name(tclSource[tc].globalOffset())
                  + "@" + rootFiles[rf]
                );
            }
            
            fileNameList subDirs
            (
                readDir(sourcePath, fileName::DIRECTORY)
            );
            
            forAll(subDirs, sd)
            {
                mkDir(destPath/subDirs[sd]);
                
                fileNameList subDirFiles
                (
                    readDir(sourcePath/subDirs[sd], fileName::FILE)
                );
                forAll(subDirFiles, sdf)
                {
                    cp
                    (
                        sourcePath/subDirs[sd]/subDirFiles[sdf],
                        destPath/subDirs[sd]/tclSource[tc].solverDomainName()
                      + "@" + name(tclSource[tc].superLoop())
                      + "@" + subDirFiles[sdf]
                      + "@" + name(tclSource[tc].globalOffset())
                    );
                } // end forAll(subDirFiles, sdf)
            } // end forAll(subDirs, sd)
        } // end forAll instants
    } // end forAll timeClusters
    if (tclSource.size())
    {
        setSolverDomainPostProcessing(tclSource[0].solverDomainName());
    }
}


void Foam::multiSolver::postCondition(const word& processor)
{
    fileName path(multiDictRegistry_.path());
    if (processor.size())
    {
        path = path/processor;
    }

    timeClusterList purgeMe(readAllTimes(processor));
    purgeMe.purgeEmpties();

    // Purge these directories
    forAll(purgeMe, i)
    {
        fileName purgePath
        (
            findInstancePath(purgeMe[i], 0).path()
        );
        rmDir(purgePath);
    }

    instantList times(Time::findTimes(path));

    forAll(times, t)
    {
        // Ignore "constant" if it exists
        if (times[t].name() == "constant")
        {
            continue;
        }
        fileName sourcePath
        (
            path/times[t].name()
        );
        
        // If timeFormat is not general, it will miss the -1 initial directory
        if (!exists(sourcePath))
        {
            if (times[t].value() == -1)
            {
                sourcePath = path/"-1";
            }
            else
            {
                continue;
            }
        }
        
        // Root files first
        fileNameList rootFiles
        (
            readDir(sourcePath, fileName::FILE)
        );
        forAll(rootFiles, rf)
        {
            timeCluster tcSubject
            (
                parseConditionedFile(rootFiles[rf], times[t])
            );

            fileName destPath;
            if (tcSubject.superLoop() == -1)
            {
                destPath = path/"multiSolver"
                    /tcSubject.solverDomainName()
                    /"initial/0";
            }
            else
            {
                destPath = path/"multiSolver"
                    /tcSubject.solverDomainName()
                    /name(tcSubject.superLoop())
                    /times[t].name();
            }
            mkDir(destPath);

            // Create multiSolverTime dictionary if it doesn't exist
            if
            (
                !exists(destPath.path()/"multiSolverTime")
             && (tcSubject.superLoop() != -1)
            )
            {
                IOdictionary multiSolverTime
                (
                    IOobject
                    (
                        "multiSolverTime",
                        multiDictRegistry_.constant(),
                        multiDictRegistry_,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE,
                        false
                    )
                );
                multiSolverTime.set("globalOffset", tcSubject.globalOffset());
                
                // Write multiSolverTime to the case/constant directory, then
                // move to destination path
                multiSolverTime.regIOobject::write();
                mv
                (
                    multiDictRegistry_.constantPath()/"multiSolverTime",
                    destPath.path()
                );
            }
            cp
            (
                sourcePath/rootFiles[rf],
                destPath/tcSubject.preConName()
            );
        } // end forAll(rootFiles, rf)

        // Subdirectories now
        fileNameList subDirs
        (
            readDir(sourcePath, fileName::DIRECTORY)
        );
        
        forAll(subDirs, sd)
        {
            fileNameList subDirFiles
            (
                readDir(sourcePath/subDirs[sd], fileName::FILE)
            );
            forAll(subDirFiles, sdf)
            {
                timeCluster tcSubject
                (
                    parseConditionedFile(subDirFiles[sdf], times[t])
                );
                fileName destPath;
                if (tcSubject.superLoop() == -1)
                {
                    destPath = sourcePath.path()/"multiSolver"
                        /tcSubject.solverDomainName()
                        /"initial/0"/subDirs[sd];
                }
                else
                {
                    destPath = sourcePath.path()/"multiSolver"
                        /tcSubject.solverDomainName()
                        /name(tcSubject.superLoop())
                        /times[t].name()
                        /subDirs[sd];
                }
                mkDir(destPath);
                cp
                (
                    sourcePath/subDirs[sd]/subDirFiles[sdf],
                    destPath/subDirs[sd]/tcSubject.preConName()
                );
            } // end forAll(subDirFiles, sdf)
        } // end forAll(subDirs, sd)
    } // end forAll(rootDirs, rd)

    // Delete root time directories
    purgeTimeDirs(path);
}


void Foam::multiSolver::setSolverDomain(const Foam::word& solverDomainName)
{
    if (run())
    {
        if (currentSolverDomain_ == "default")
        {
            setInitialSolverDomain(solverDomainName);
        }
        else
        {
            setNextSolverDomain(solverDomainName);
        }
    }
    if (Pstream::parRun())
    {
        synchronizeParallel();
    }
}


void Foam::multiSolver::setSolverDomainPostProcessing
(
    const Foam::word& solverDomainName
)
{
    if (!solverDomains_.found(solverDomainName))
    {
        FatalErrorIn("multiSolver::setInitialSolverDomain")
            << "Initial solverDomainName '" << solverDomainName << "' does"
            << " not exist in multiSolver dictionary.  Found entries are: "
            << solverDomains_.toc()
            << abort(FatalError);
    }

    currentSolverDomain_ = solverDomainName;
    
    setSolverDomainControls(currentSolverDomain_);

    // startTime is set to the earliest in case/[time]
    instantList il(multiDictRegistry_.times());
    label first(Time::findClosestTimeIndex(il,-1.0));

    startTime_ = il[first].value();

    word stopAtSetting("endTime");

    // Build the new controlDict
    IOdictionary newControlDict
    (
        IOobject
        (
            Time::controlDictName,
            multiDictRegistry_.system(),
            multiDictRegistry_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE,
            false
        ),
        currentSolverDomainDict_
    );
    
    // Remove multiSolver-specific values from dictionary
    newControlDict.remove("startFrom");
    newControlDict.remove("startTime");
    newControlDict.remove("stopAt");
    newControlDict.remove("endTime");
    newControlDict.remove("iterations");
    newControlDict.remove("purgeWriteSuperLoops");
    newControlDict.remove("timeFormat");
    newControlDict.remove("timePrecision");
    newControlDict.remove("storeFields");
    newControlDict.remove("elapsedTime");
    
    // Add values to obtain the desired behaviour
    newControlDict.set("startFrom", "startTime");
    newControlDict.set("startTime", startTime_);
    newControlDict.set("stopAt", stopAtSetting);
    newControlDict.set("endTime", endTime_);
    if (multiSolverControl_.found("timeFormat"))
    {    
        newControlDict.set
        (
            "timeFormat",
            word(multiSolverControl_.lookup("timeFormat"))
        );
    }
    if (multiSolverControl_.found("timePrecision"))
    {    
        newControlDict.set
        (
            "timePrecision",
            readScalar(multiSolverControl_.lookup("timePrecision"))
        );
    }

    // Write the dictionary to the case directory
    newControlDict.regIOobject::write();

    swapDictionaries(currentSolverDomain_);
}


Foam::multiSolver& Foam::multiSolver::operator++()
{
    superLoop_++;
    noSaveSinceSuperLoopIncrement_ = true;
    return *this;
}


Foam::multiSolver& Foam::multiSolver::operator++(int)
{
    return operator++();
}


bool Foam::multiSolver::run() const
{
    // If case/[time] are present, run must continue to next 'setSolverDomain'
    // so that they are archived properly.
    instantList il(Time::findTimes(multiDictRegistry_.path()));
    return !(finished_ && (il.size() == 1));
}


bool Foam::multiSolver::end() const
{
    // If case/[time] are present, run must continue to next 'setSolverDomain'
    // so that they are archived properly.
    instantList il(Time::findTimes(multiDictRegistry_.path()));
    return (finished_ && (il.size() == 1));
}

#include "multiSolverSetControls.C"
#include "multiSolverSetInitialSolverDomain.C"
#include "multiSolverSetNextSolverDomain.C"
#include "multiSolverTimeFunctions.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

