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

#include "argList.H"
#include "OSspecific.H"
#include "clock.H"
#include "IFstream.H"
#include "dictionary.H"
#include "Switch.H"
#include "IOobject.H"
#include "JobInfo.H"
#include "labelList.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    SLList<string>    argList::validArgs;
    HashTable<string> argList::validOptions;
    HashTable<string> argList::validParOptions;
}


Foam::argList::initValidTables::initValidTables()
{
    validOptions.insert("case", "dir");
    validOptions.insert("parallel", "");
    validParOptions.insert("parallel", "");

    Pstream::addValidParOptions(validParOptions);
}


Foam::argList::initValidTables dummyInitValidTables;


// convert argv -> args_
// transform sequences with "(" ... ")" into string lists in the process
bool Foam::argList::regroupArgv(int& argc, char**& argv)
{
    int level = 0;
    int nArgs = 0;
    string tmpString;

    // note: we also re-write directly into args_
    // and use a second pass to sort out args/options
    for (int argi=0; argi < argc; argi++)
    {
        if (strcmp(argv[argi], "(") == 0)
        {
            level++;
            tmpString += "(";
        }
        else if (strcmp(argv[argi], ")") == 0)
        {
            if (level >= 1)
            {
                level--;
                tmpString += ")";
                if (level == 0)
                {
                    args_[nArgs++] = tmpString;
                    tmpString.clear();
                }
            }
            else
            {
                args_[nArgs++] = argv[argi];
            }
        }
        else if (level)
        {
            // quote each string element
            tmpString += "\"";
            tmpString += argv[argi];
            tmpString += "\"";
        }
        else
        {
            args_[nArgs++] = argv[argi];
        }
    }

    if (tmpString.size())
    {
        args_[nArgs++] = tmpString;
    }

    args_.setSize(nArgs);

    return nArgs < argc;
}


// get rootPath_ / globalCase_ from one of the following forms
//   * [-case dir]
//   * cwd
void Foam::argList::getRootCase()
{
    fileName casePath;

    // [-case dir] specified
    HashTable<string>::iterator iter = options_.find("case");

    if (iter != options_.end())
    {
        casePath = iter();
        casePath.removeRepeated('/');
        casePath.removeTrailing('/');
    }
    else
    {
        // nothing specified, use the current dir
        casePath = cwd();

        // we could add this back in as '-case'?
        // options_.insert("case", casePath);
    }

    rootPath_   = casePath.path();
    globalCase_ = casePath.name();
    case_       = globalCase_;
}


Foam::stringList::subList Foam::argList::additionalArgs() const
{
    return stringList::subList(args_, args_.size() - 1, 1);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::argList::argList
(
    int& argc,
    char**& argv,
    bool checkArgs,
    bool checkOpts
)
:
    args_(argc),
    options_(argc)
{
    // Check if this run is a parallel run by searching for any parallel option
    // If found call runPar (might filter argv)

    for (int argi=0; argi<argc; argi++)
    {
        if (argv[argi][0] == '-')
        {
            const char *optionName = &argv[argi][1];

            if (validParOptions.found(optionName))
            {
                parRunControl_.runPar(argc, argv);
                break;
            }
        }
    }

    // Print the banner once only for parallel runs
    if (Pstream::master())
    {
        IOobject::writeBanner(Info, true);
    }

    // convert argv -> args_ and capture ( ... ) lists
    // for normal arguments and for options
    regroupArgv(argc, argv);

    // Get executable name
    args_[0]    = fileName(argv[0]);
    executable_ = fileName(argv[0]).name();

    // Check arguments and options, we already have argv[0]
    int nArgs = 1;
    string argListString = args_[0];

    for (int argi = 1; argi < args_.size(); argi++)
    {
        argListString += ' ';
        argListString += args_[argi];

        if (args_[argi][0] == '-')
        {
            const char *optionName = &args_[argi][1];

            if
            (
                (
                    validOptions.found(optionName)
                 && validOptions[optionName] != ""
                )
             || (
                    validParOptions.found(optionName)
                 && validParOptions[optionName] != ""
                )
            )
            {
                argi++;
                if (argi >= args_.size())
                {
                    FatalError
                        << "option " << "'-" << optionName << '\''
                        << " requires an argument"
                        << exit(FatalError);
                }

                argListString += ' ';
                argListString += args_[argi];
                options_.insert(optionName, args_[argi]);
            }
            else
            {
                options_.insert(optionName, "");
            }
        }
        else
        {
            if (nArgs != argi)
            {
                args_[nArgs] = args_[argi];
            }
            nArgs++;
        }
    }

    args_.setSize(nArgs);

    // Help/documentation options:
    //   -help    print the usage
    //   -doc     display application documentation in browser
    //   -srcDoc  display source code in browser
    if
    (
        options_.found("help")
     || options_.found("doc")
     || options_.found("srcDoc")
    )
    {
        if (options_.found("help"))
        {
            printUsage();
        }

        // only display one or the other
        if (options_.found("srcDoc"))
        {
            displayDoc(true);
        }
        else if (options_.found("doc"))
        {
            displayDoc(false);
        }

        ::exit(0);
    }

    // Print the usage message and exit if the number of arguments is incorrect
    if (!check(checkArgs, checkOpts))
    {
        FatalError.exit();
    }


    string dateString = clock::date();
    string timeString = clock::clockTime();

    if (Pstream::master())
    {
        Info<< "Exec   : " << argListString.c_str() << nl
            << "Date   : " << dateString.c_str() << nl
            << "Time   : " << timeString.c_str() << nl
            << "Host   : " << hostName() << nl
            << "PID    : " << pid() << endl;
    }

    jobInfo.add("startDate", dateString);
    jobInfo.add("startTime", timeString);
    jobInfo.add("userName", userName());
    jobInfo.add("foamVersion", word(FOAMversion));
    jobInfo.add("code", executable_);
    jobInfo.add("argList", argListString);
    jobInfo.add("currentDir", cwd());
    jobInfo.add("PPID", ppid());
    jobInfo.add("PGID", pgid());


    // Case is a single processor run unless it is running parallel
    int nProcs = 1;

    // If this actually is a parallel run
    if (parRunControl_.parRun())
    {
        // For the master
        if (Pstream::master())
        {
            fileNameList roots;

            // establish rootPath_/globalCase_/case_ for master
            getRootCase();

            IFstream decompDictStream
            (
                rootPath_/globalCase_/"system/decomposeParDict"
            );

            if (!decompDictStream.good())
            {
                FatalError
                    << "Cannot read "
                    << decompDictStream.name()
                    << exit(FatalError);
            }

            dictionary decompositionDict(decompDictStream);

            Switch distributed(false);

            if (decompositionDict.found("distributed"))
            {
                decompositionDict.lookup("distributed") >> distributed;

                if (distributed)
                {
                    decompositionDict.lookup("roots") >> roots;

                    if (roots.size() != Pstream::nProcs())
                    {
                        FatalError
                            << "number of entries in "
                            << "decompositionDict::roots"
                            << " is not equal to the number of processors "
                            << Pstream::nProcs()
                            << exit(FatalError);
                    }
                }
            }


            label dictNProcs
            (
                readLabel
                (
                    decompositionDict.lookup("numberOfSubdomains")
                )
            );

            // Check number of processors. We have nProcs(number of
            // actual processes), dictNProcs(wanted number of processes read
            // from decompositionDict) and nProcDirs(number of processor
            // directories - n/a when running distributed)
            //
            // - normal running : nProcs = dictNProcs = nProcDirs
            // - decomposition to more processors : nProcs = dictNProcs
            // - decomposition to less processors : nProcs = nProcDirs
            if (dictNProcs > Pstream::nProcs())
            {
                FatalError
                    << decompDictStream.name()
                    << " specifies " << dictNProcs
                    << " processors but job was started with "
                    << Pstream::nProcs() << " processors."
                    << exit(FatalError);
            }

            if (!distributed && dictNProcs < Pstream::nProcs())
            {
                // Possibly going to fewer processors.
                // Check if all procDirs are there.
                label nProcDirs = 0;
                while
                (
                    dir
                    (
                        rootPath_/globalCase_/"processor"
                      + name(++nProcDirs)
                    )
                )
                {}

                if (nProcDirs != Pstream::nProcs())
                {
                    FatalError
                        << "number of processor directories = "
                        << nProcDirs
                        << " is not equal to the number of processors = "
                        << Pstream::nProcs()
                        << exit(FatalError);
                }
            }

            // distributed data
            if (roots.size())
            {
                bool hadOptCase = options_.found("case");

                // Distribute the master's argument list (with new root)
                for
                (
                    int slave=Pstream::firstSlave();
                    slave<=Pstream::lastSlave();
                    slave++
                )
                {
                    options_.erase("case");
                    options_.insert
                    (
                        "case",
                        fileName(roots[slave-1])/globalCase_
                    );

                    OPstream toSlave(Pstream::scheduled, slave);
                    toSlave << args_ << options_;
                }

                options_.erase("case");

                // restore [-case dir]
                if (hadOptCase)
                {
                    options_.insert("case", rootPath_/globalCase_);
                }
            }
            else
            {
                // Distribute the master's argument list (unaltered)
                for
                (
                    int slave=Pstream::firstSlave();
                    slave<=Pstream::lastSlave();
                    slave++
                )
                {
                    OPstream toSlave(Pstream::scheduled, slave);
                    toSlave << args_ << options_;
                }
            }
        }
        else
        {
            // Collect the master's argument list
            IPstream fromMaster(Pstream::scheduled, Pstream::masterNo());
            fromMaster >> args_ >> options_;

            // establish rootPath_/globalCase_/case_ for slave
            getRootCase();
        }

        nProcs = Pstream::nProcs();
        case_ = globalCase_/(word("processor") + name(Pstream::myProcNo()));
    }
    else
    {
        // establish rootPath_/globalCase_/case_
        getRootCase();

        case_ = globalCase_;
    }


    wordList slaveProcs;

    // collect slave machine/pid
    if (parRunControl_.parRun())
    {
        if (Pstream::master())
        {
            slaveProcs.setSize(Pstream::nProcs() - 1);
            word  slaveMachine;
            label slavePid;

            label procI = 0;
            for
            (
                int slave=Pstream::firstSlave();
                slave<=Pstream::lastSlave();
                slave++
            )
            {
                IPstream fromSlave(Pstream::scheduled, slave);
                fromSlave >> slaveMachine >> slavePid;

                slaveProcs[procI++] = slaveMachine + "." + name(slavePid);
            }
        }
        else
        {
            OPstream toMaster(Pstream::scheduled, Pstream::masterNo());
            toMaster << hostName() << pid();
        }
    }


    if (Pstream::master())
    {
        Info<< "Case   : " << (rootPath_/globalCase_).c_str() << nl
            << "nProcs : " << nProcs << endl;
    }

    if (parRunControl_.parRun() && Pstream::master())
    {
        Info<< "Slaves : " << slaveProcs << nl
            << "Pstream initialized with:" << nl
            << "    floatTransfer     : " << Pstream::floatTransfer << nl
            << "    nProcsSimpleSum   : " << Pstream::nProcsSimpleSum << nl
            << "    commsType         : "
            << Pstream::commsTypeNames[Pstream::defaultCommsType]
            << endl;
    }

    jobInfo.add("root", rootPath_);
    jobInfo.add("case", globalCase_);
    jobInfo.add("nProcs", nProcs);
    if (slaveProcs.size())
    {
        jobInfo.add("slaves", slaveProcs);
    }
    jobInfo.write();

    // Set the case as an environment variable
    setEnv("FOAM_CASE", rootPath_/globalCase_, true);

    // Switch on signal trapping. We have to wait until after Pstream::init
    // since this sets up its own ones.
    sigFpe_.set();
    sigInt_.set();
    sigQuit_.set();
    sigSegv_.set();

    if (Pstream::master())
    {
        Info<< endl;
        IOobject::writeDivider(Info);
    }
}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

Foam::argList::~argList()
{
    jobInfo.end();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::argList::noParallel()
{
    validOptions.erase("parallel");
}


void Foam::argList::printUsage() const
{
    Info<< nl
        << "Usage: " << executable_;

    for
    (
        SLList<string>::iterator iter = validArgs.begin();
        iter != validArgs.end();
        ++iter
    )
    {
        Info<< " <" << iter().c_str() << '>';
    }

    for
    (
        HashTable<string>::iterator iter = validOptions.begin();
        iter != validOptions.end();
        ++iter
    )
    {
        Info<< " [-" << iter.key();

        if (iter().size())
        {
            Info<< ' ' << iter().c_str();
        }

        Info<< ']';
    }

    // place help/doc options of the way at the end,
    // but with an extra space to separate it a little
    Info<< "  [-help] [-doc] [-srcDoc]" << endl;
}


void Foam::argList::displayDoc(bool source) const
{
    const dictionary& docDict = debug::controlDict().subDict("Documentation");
    List<fileName> docDirs(docDict.lookup("doxyDocDirs"));
    List<fileName> docExts(docDict.lookup("doxySourceFileExts"));

    // for source code: change foo_8C.html to foo_8C-source.html
    if (source)
    {
        forAll(docExts, extI)
        {
            docExts[extI].replace(".", "-source.");
        }
    }

    fileName docFile;
    bool found = false;

    forAll(docDirs, dirI)
    {
        forAll(docExts, extI)
        {
            docFile = docDirs[dirI]/executable_ + docExts[extI];
            docFile.expand();

            if (exists(docFile))
            {
                found = true;
                break;
            }
        }
        if (found)
        {
            break;
        }
    }

    if (found)
    {
        string docBrowser(docDict.lookup("docBrowser"));
        docBrowser.replaceAll("%f", docFile);

        Info<< "Show documentation: " << docBrowser.c_str() << endl;

        system(docBrowser);
    }
    else
    {
        Info<< "No documentation found" << endl;
    }
}


bool Foam::argList::check(bool checkArgs, bool checkOpts) const
{
    bool ok = true;

    if (Pstream::master())
    {
        if (checkArgs && args_.size() - 1 != validArgs.size())
        {
            FatalError
                << "Wrong number of arguments, expected " << validArgs.size()
                << " found " << args_.size() - 1 << endl;
            ok = false;
        }

        if (checkOpts)
        {
            forAllConstIter(HashTable<string>, options_, iter)
            {
                if
                (
                    !validOptions.found(iter.key())
                 && !validParOptions.found(iter.key())
                )
                {
                    FatalError
                        << "Invalid option: -" << iter.key() << endl;
                    ok = false;
                }
            }
        }

        if (!ok)
        {
            printUsage();
        }
    }

    return ok;
}


bool Foam::argList::checkRootCase() const
{
    if (!dir(rootPath()))
    {
        FatalError
            << executable_
            << ": cannot open root directory " << rootPath()
            << endl;

        return false;
    }

    if (!dir(path()) && Pstream::master())
    {
        // Allow slaves on non-existing processor directories, created later
        FatalError
            << executable_
            << ": cannot open case directory " << path()
            << endl;

        return false;
    }

    return true;
}


// ************************************************************************* //
