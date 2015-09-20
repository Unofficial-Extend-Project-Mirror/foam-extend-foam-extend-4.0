/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Application
    multiSolver

Description
    Post-processor utility required for working with multiSolver-enabled
    applications.  These applications store data output in a different
    location than usual.  This utility loads and unloads data from that
    location to where post-processors expect it.

Author
    David L. F. Gaden

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "multiSolver.H"

void parseOptions
(
    wordList * ptrSolverDomains,
    labelList * ptrSuperLoops,
    string options
)
{
    IStringStream optionsStream(options);
    label nSolverDomains(0);
    label nSuperLoops(0);

    // Get solverDomainNames, if any
    while (not optionsStream.eof())
    {
        token nextOption(optionsStream);

        // Bug workaround
        if (nextOption.type() == token::FATALERROR)
        {
            break;
        }

        if (nextOption.isLabel())
        {
            ptrSuperLoops->setSize(++nSuperLoops);
            ptrSuperLoops->operator[](nSuperLoops - 1)
                = nextOption.labelToken();
            break;
        }
        if (nextOption.isWord())
        {
            ptrSolverDomains->setSize(++nSolverDomains);
            ptrSolverDomains->operator[](nSolverDomains - 1)
                = nextOption.wordToken();
        }
        else
        {
            // not word, not label, fail
            FatalErrorIn("multiSolver::parseOptions")
                << "Expecting word or label.  Neither found at position "
                << nSolverDomains - 1 << " in " << options
                << abort(FatalError);
        }
    }

    // Get superLoopList
    while (not optionsStream.eof())
    {
        token nextOption(optionsStream);

        // Bug workaround
        if (nextOption.type() == token::FATALERROR)
        {
            break;
        }

        if (nextOption.isLabel())
        {
            ptrSuperLoops->setSize(++nSuperLoops);
            ptrSuperLoops->operator[](nSuperLoops - 1)
                = nextOption.labelToken();
        }
        else if (nSuperLoops > 0)
        {
            // might be a range -> label : label

            if (nextOption.isPunctuation())
            {
                token::punctuationToken p(nextOption.pToken());
                if (p == token::COLON)
                {
                    token nextNextOption(optionsStream);
                    if (nextNextOption.isLabel())
                    {
                        label toValue(nextNextOption.labelToken());
                        label fromValue
                        (
                            ptrSuperLoops->operator[](nSuperLoops - 1)
                        );

                        if (toValue > fromValue)
                        {
                            // correct range format
                            for (label i = fromValue + 1; i <= toValue; i++)
                            {
                                ptrSuperLoops->setSize(++nSuperLoops);
                                ptrSuperLoops->operator[](nSuperLoops - 1) = i;
                            }
                        }
                        else
                        {
                            // greater than / less than, range
                            FatalErrorIn("multiSolver::parseOptions")
                                << "superLoop range incorrect order.  'from : "
                                << "to' where 'from' should be less than "
                                << "'to'.  Values read are '" << fromValue
                                << " : " << toValue
                                << abort(FatalError);
                        }
                    }
                    else
                    {
                        // nextNext not label
                        FatalErrorIn("multiSolver::parseOptions")
                            << "Incorrect syntax.  Expecting label after ':' "
                            << "in " << options
                            << abort(FatalError);
                    }
                }
                else
                {
                    // non : punctuation
                    FatalErrorIn("multiSolver::parseOptions")
                        << "Incorrect syntax.  Expecting label, word, or ':' "
                        << "in " << options
                        << abort(FatalError);
                }
            }
            else
            {
                // not punctuation
                FatalErrorIn("multiSolver::parseOptions")
                    << "Incorrect syntax.  Expecting label, word, or ':' "
                    << "in " << options
                    << abort(FatalError);
            }
        }
        else
        {
            // not label, not word
            FatalErrorIn("multiSolver::parseOptions")
                << "Incorrect syntax.  Expecting label, word, or ':' "
                << "in " << options
                << abort(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::validOptions.insert("list","");
    argList::validOptions.insert
    (
        "load","<[solverDomainName] [superLoopNumber(s)]>"
    );
    argList::validOptions.insert
    (
        "purge","<[solverDomainName] [superLoopNumber(s)]>"
    );
    argList::validOptions.insert("set","<solverDomainName>");
    argList::validOptions.insert("preDecompose", "");
    argList::validOptions.insert("postDecompose", "");
    argList::validOptions.insert("preReconstruct", "");
    argList::validOptions.insert("postReconstruct", "");

    argList::validOptions.insert("global","");
    argList::validOptions.insert("local","");

    // default behaviour is purge the case/[time] directory before '-load'
    // command.  '-noPurge' prevents this.  Allows for more complicated
    // load data selections by executing multiSolver several times
    argList::validOptions.insert("noPurge","");

    // default behaviour is: if there is only one solverDomain specified, use
    // setSolverDomain() on it.  Same as multiSolver -set solverDomain.
    // '-noSet' prevents this.
    argList::validOptions.insert("noSet","");

    // default behaviour is: if there are storeFields defined, when loading, it
    // will copy the store fields into every time instance where they are
    // absent.  '-noStore' will prevent this.
    argList::validOptions.insert("noStore","");

#   include "setRootCase.H"

    enum commandType
    {
        list,
        load,
        purge,
        set,
        preDecompose,
        postDecompose,
        preReconstruct,
        postReconstruct
    };
    commandType command;
    string options;
    bool global = false;
    bool local = false;
    bool all = false;
    bool root = false;
    bool noPurge = false;
    bool noSet = false;
    bool noStore = false;
    label nCommands(0);

    // Read arguments
    if (args.optionFound("list"))
    {
        nCommands++;
        command = list;
    }
    if (args.optionFound("load"))
    {
        nCommands++;
        command = load;
        options = args.options()["load"];
    }
    if (args.optionFound("purge"))
    {
        nCommands++;
        command = purge;
        options = args.options()["purge"];
    }
    if (args.optionFound("set"))
    {
        nCommands++;
        command = set;
        options = args.options()["set"];
    }
    if (args.optionFound("preDecompose"))
    {
        nCommands++;
        command = preDecompose;
    }
    if (args.optionFound("postDecompose"))
    {
        nCommands++;
        command = postDecompose;
    }
    if (args.optionFound("preReconstruct"))
    {
        nCommands++;
        command = preReconstruct;
    }
    if (args.optionFound("postReconstruct"))
    {
        nCommands++;
        command = postReconstruct;
    }
    if (args.optionFound("global"))
    {
        global = true;
    }
    if (args.optionFound("local"))
    {
        local = true;
    }
    if (args.optionFound("noPurge"))
    {
        noPurge = true;
    }
    if (args.optionFound("noSet"))
    {
        noSet = true;
    }
    if (args.optionFound("noStore"))
    {
        noStore = true;
    }

    // Error checking
    if (nCommands == 0)
    {
        FatalErrorIn("multiSolver::main")
            << "multiSolver - nothing to do.  Use 'multiSolver -help' for assistance."
            << abort(FatalError);
    }
    else if (nCommands > 1)
    {
        FatalErrorIn("multiSolver::main")
            << "More than one command found.  Use only one of:"
            << "\n\t-list"
            << "\n\t-purge"
            << "\n\t-set"
            << "\n\t-preDecompose"
            << "\n\t-postDecompose"
            << "\n\t-preReconstruct"
            << "\n\t-postReconstruct\n"
            << abort(FatalError);
    }
    if (global && local)
    {
        FatalErrorIn("multiSolver::main")
            << "Options global and local both specified.  Use only one or "
            << "none."
            << abort(FatalError);
    }
    if ((command != load) && (noPurge || noSet || noStore))
    {
        FatalErrorIn("multiSolver::main")
            << "'noPurge', 'noSet' and 'noStore' can only be used with the "
            << "'-load' command."
            << abort(FatalError);
    }

    multiSolver multiRun
    (
        Foam::multiSolver::multiControlDictName,
        args.rootPath(),
        args.caseName()
    );

    const IOdictionary& mcd(multiRun.multiControlDict());
    wordList solverDomains(0);
    labelList superLoops(0);
    if
    (
        (command != list)
     && (command != preDecompose)
     && (command != postDecompose)
     && (command != preReconstruct)
     && (command != postReconstruct)
    )
    {
        parseOptions(&solverDomains, &superLoops, options);
    }

    // Special words - all, root
    if (solverDomains.size() == 1)
    {
        if (solverDomains[0] == "all")
        {
            all = true;
        }
        else if (solverDomains[0] == "root")
        {
            root = true;
        }
    }

    // More error checking
    if (root && ((command == load) || (command == set)))
    {
        FatalErrorIn("multiSolver::main")
            << "'root' is not a valid option with '-load' or '-set'"
            << abort(FatalError);
    }
    if (all && (command == set))
    {
        FatalErrorIn("multiSolver::main")
            << "'all' is not a valid option with '-set'"
            << abort(FatalError);
    }
    if ((command == set) && ((solverDomains.size() > 1) || superLoops.size()))
    {
        FatalErrorIn("multiSolver::main")
            << "'-set' can only have a single solverDomain name as an option."
            << abort(FatalError);
    }
    if (all && superLoops.size())
    {
        FatalErrorIn("multiSolver::main")
            << "'all' cannot be followed by superLoop numbers.  To specify "
            << "a superLoop range for all solverDomains, omit the solverDomain"
            << " name entirely.  e.g. multiSolver -load '0:4 6'"
            << abort(FatalError);
    }
    if (root && superLoops.size())
    {
        FatalErrorIn("multiSolver::main")
            << "'root' cannot be followed by superLoop numbers.  'root' refers"
            << " to case/[time] directories.  There are no superLoops here."
            << abort(FatalError);
    }

    // Check for correct solverDomain names
    if
    (
        !all
     && !root
     && (command != list)
     && (command != preDecompose)
     && (command != postDecompose)
     && (command != preReconstruct)
     && (command != postReconstruct)
    )
    {
        forAll(solverDomains, i)
        {
            if (solverDomains[i] == "default")
            {
                // default not permitted
                FatalErrorIn("multiSolver::main")
                    << "'default' is not a permitted solverDomain name."
                    << abort(FatalError);
            }
            if (!mcd.subDict("solverDomains").found(solverDomains[i]))
            {
                // Incorrect solver domain name
                FatalErrorIn("multiSolver::main")
                    << "solverDomainName " << solverDomains[i] << "is not "
                    << "found."
                    << abort(FatalError);
            }
        }
    }

    // Load specified timeClusterLists
    timeClusterList tclSource(0);

    if (all)
    {
        // read all
        tclSource = multiRun.readAllTimes();
        forAll(tclSource, i)
        {
            if (tclSource[i].superLoop() == -1)
            {
                tclSource[i].times().clear();
            }
        }
        tclSource.purgeEmpties();
    }
    else if
    (
        !superLoops.size()
     && (command != set)
     && (command != list)
     && (command != preDecompose)
     && (command != postDecompose)
     && (command != preReconstruct)
     && (command != postReconstruct)
    )
    {
        // no superLoops specified - read entire solverDomains
        forAll (solverDomains, sd)
        {
            tclSource.append
            (
                multiRun.readSolverDomainTimes(solverDomains[sd])
            );
        }
    }
    else if
    (
        !root
     && (command != set)
     && (command != list)
     && (command != preDecompose)
     && (command != postDecompose)
     && (command != preReconstruct)
     && (command != postReconstruct)
    )
    {
        // read individual superLoops
        if (!solverDomains.size())
        {
            solverDomains = mcd.subDict("solverDomains").toc();
        }
        forAll(superLoops, sl)
        {
            forAll(solverDomains, sd)
            {
                if (solverDomains[sd] == "default") continue;
                tclSource.append
                (
                    multiRun.readSuperLoopTimes
                    (
                        solverDomains[sd],
                        superLoops[sl]
                    )
                );
            }
        }
    }

    if (tclSource.size())
    {
        if (!tclSource.purgeEmpties())
        {
            FatalErrorIn("multiSolver::main")
                << "No data found with specified parameters."
                << abort(FatalError);
        }
    }

    switch (command)
    {
        case list:
        {
            Info << "Listing available data:\n" << endl;
            Info << "superLoops by solverDomain:" << endl;
            solverDomains = mcd.subDict("solverDomains").toc();
            fileName listPath
            (
                multiRun.multiDictRegistry().path()/"multiSolver"
            );

            forAll(solverDomains, i)
            {
                if (solverDomains[i] == "default") continue;
                Info << solverDomains[i] << ":" << endl;
                Info << multiRun.findSuperLoops(listPath/solverDomains[i])
                    << endl;
            }
            Info << endl;
            break;
        }
        case load:
        {
            // Default behaviour - use local time unless overlapping, then use
            // global time; if global overlaps, fail.  -local and -global force
            // the behaviour
            bool localOverlap(!multiRun.nonOverlapping(tclSource, false));
            bool globalOverlap(!multiRun.nonOverlapping(tclSource, true));
            if (local && localOverlap)
            {
                FatalErrorIn("multiSolver::main")
                    << "'-local' option used for data with overlapping local "
                    << "values.  Try using a single solverDomain / superLoop, "
                    << "or leave '-local' off."
                    << abort(FatalError);
            }
            if (globalOverlap)
            {
                FatalErrorIn("multiSolver::main")
                    << "globalTime values are overlapping.  This should not "
                    << "happen.  Ensure you have not specified the same "
                    << "solverDomain and/or superLoop more than once.  If "
                    << "that fails, try using 'multiSolver -purge all' and "
                    << "rerunning the simulation.  If the problem persists, "
                    << "it is a bug."
                    << abort(FatalError);
            }

            if (!noPurge)
            {
                Info << "Purging existing time directories in case root"
                    << endl;
                multiRun.purgeTimeDirs(multiRun.multiDictRegistry().path());
            }

            Info << "Loading data from multiSolver directories to case root"
                << endl;
            if
            (
                !multiRun.loadTimeClusterList
                (
                    tclSource,
                    global || localOverlap,
                    !noStore
                )
            )
            {
                FatalErrorIn("multiRun::main")
                    << "loadTimeClusterList failed.  timeClusterList contents: "
                    << tclSource
                    << abort(FatalError);
            }
            break;
        }
        case purge:
            if (root)
            {
                Info << "Purging time directories from case root" << endl;
                multiRun.purgeTimeDirs(multiRun.multiDictRegistry().path());
            }
            else
            {
                Info << "Purging time directories from multiSolver directories"
                    << endl;
                forAll(tclSource, i)
                {
                    // do not purge 'initial' directory, even if specified
                    if (tclSource[i].superLoop() < 0) continue;
                    fileName purgePath
                    (
                        multiRun.findInstancePath(tclSource[i], 0).path()
                    );
                    rmDir(purgePath);
                }
            }
            break;
        case set:
            // do nothing here
            break;
        case preDecompose:
        {
            Info << "Performing preDecompose" << endl;
            multiRun.preCondition();
            break;
        }
        case postDecompose:
        {
            Info << "Performing postDecompose" << endl;

            fileNameList dirEntries
            (
                readDir
                (
                    multiRun.multiDictRegistry().path(), fileName::DIRECTORY
                )
            );

            forAll(dirEntries, de)
            {
                if (dirEntries[de](9) == "processor")
                {
                    Info << "Reading " << dirEntries[de] << endl;

                    multiRun.postCondition
                    (
                        dirEntries[de]
                    );

                    // Copy system to processorN
                    cp
                    (
                        multiRun.multiDictRegistry().path()
                            /multiRun.multiDictRegistry().system(),
                        multiRun.multiDictRegistry().path()/dirEntries[de]
                    );

                    // Copy constant/files to processorN/constant
                    fileNameList constantContents
                    (
                        readDir
                        (
                            multiRun.multiDictRegistry().path()
                                /multiRun.multiDictRegistry().constant(),
                            fileName::FILE
                        )
                    );
                    forAll(constantContents, cc)
                    {
                        cp
                        (
                            multiRun.multiDictRegistry().path()
                                /multiRun.multiDictRegistry().constant()
                                /constantContents[cc],
                            multiRun.multiDictRegistry().path()/dirEntries[de]
                                /multiRun.multiDictRegistry().constant()
                        );
                    }

                    // Copy constant/directories to processorN/constant
                    constantContents = readDir
                    (
                        multiRun.multiDictRegistry().path()
                            /multiRun.multiDictRegistry().constant(),
                        fileName::DIRECTORY
                    );
                    forAll(constantContents, cc)
                    {
                        // Ingore mesh directory
                        if (constantContents[cc] == "polyMesh")
                        {
                            continue;
                        }
                        cp
                        (
                            multiRun.multiDictRegistry().path()
                                /multiRun.multiDictRegistry().constant()
                                /constantContents[cc],
                            multiRun.multiDictRegistry().path()/dirEntries[de]
                                /multiRun.multiDictRegistry().constant()
                        );
                    }
                }
            }
            multiRun.purgeTimeDirs(multiRun.multiDictRegistry().path());
            break;
        }
        case preReconstruct:
        {
            Info << "Performing preReconstruct" << endl;
            fileNameList dirEntries
            (
                readDir
                (
                    multiRun.multiDictRegistry().path(), fileName::DIRECTORY
                )
            );

            forAll(dirEntries, de)
            {
                if (dirEntries[de](9) == "processor")
                {
                    Info << "Reading " << dirEntries[de] << endl;
                    multiRun.preCondition
                    (
                        dirEntries[de]
                    );
                    // Fix missing 0.00000e+00 directory if it exists
                    mkDir
                    (
                        multiRun.multiDictRegistry().path()/dirEntries[de]/"0"
                    );
                }
            }
            break;
        }
        case postReconstruct:
        {
            Info << "Performing postReconstruct" << endl;

            Info << "Reading preconditioned time directories" << endl;
            multiRun.postCondition();

            Info << "Purging preconditioned time directories"
                << endl;

            // Clean up extra time directories
            fileNameList dirEntries
            (
                readDir
                (
                    multiRun.multiDictRegistry().path(), fileName::DIRECTORY
                )
            );

            forAll(dirEntries, de)
            {
                if (dirEntries[de](9) == "processor")
                {
                    multiRun.purgeTimeDirs
                    (
                        multiRun.multiDictRegistry().path()/dirEntries[de]
                    );
                }
            }
            break;
        }
    }

    // Execute set command - either from an explicit '-set' or from a '-load'
    // with only one solverDomain as an option

    if
    (
        (command == set)
     || (
            (command == load)
         && (solverDomains.size() == 1)
         && (!all)
        )
    )
    {
        Info << "Changing to " << solverDomains[0] << " settings." << endl;
        multiRun.setSolverDomainPostProcessing(solverDomains[0]);
    }

    Info << "\nCommand completed successfully.\n" << endl;
    return(0);
}

// ************************************************************************* //
