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

\*---------------------------------------------------------------------------*/


void Foam::multiSolver::setInitialSolverDomain(const word& solverDomainName)
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

    // Purge all time directories from case directory root
    purgeTimeDirs(multiDictRegistry_.path());

    // Purge any constant/superLoopData/
    fileName superLoopDataPath
    (
        multiDictRegistry_.path()/multiDictRegistry_.constant()
            /"superLoopData"
    );
    if (exists(superLoopDataPath))
    {
        rmDir(superLoopDataPath);
    }

    // Read initial settings and determine data source (from which path the
    // initial data is copied, the starting superLoop_, and the current
    // globalTime (used to determine globalOffset).  Rules that are applied:
    //
    //   1. superLoop_ = data source superLoop
    //       a. unless data source solverDomain != currentSolverDomain_, in
    //          which case, superLoop_ = data source superLoop + 1
    //   2. globalTime = data source globalTime.  globalTime does not increment
    //      when swapping solver domains.
    //   3. startTime = data source local time
    //       a. unless data source solverDomain != currentSolverDomain_, in
    //          which case, startTime is dictated by the solverDomains
    //          subdictionary.
    //   4. endTime is determined by the solverDomains subdictionary
    //       a. unless the finalStopAt trumps it

    // Find initial data source
    timeCluster tcSource(initialDataSource());

    fileName sourcePath(findInstancePath(tcSource, tcSource.size() - 1));
    superLoop_ = tcSource.superLoop();
    globalIndex_ = tcSource.globalIndex();

    // If starting from initial conditions, superLoop_ = -1
    if (superLoop_ < 0) superLoop_ = 0;
    scalar globalTime(tcSource.globalValue(tcSource.size() - 1));
    scalar localStartTime(tcSource.localValue(tcSource.size() -1));

    // Now to apply the exceptions if currentSolverDomain_ != data source
    // solverDomain (see long comment above).
    if (sourcePath.path().path().name() != currentSolverDomain_)
    {
        superLoop_++;
        globalIndex_++;

        switch (startFrom_)
        {
            case mtsFirstTime:
                localStartTime = 0;
                break;
            case mtsStartTime:
                localStartTime = startTime_;
                break;
            case mtsLatestTimeThisDomain:
                {
                    timeCluster tcTemp
                    (
                        findLatestLocalTime
                        (
                            readSolverDomainTimes(currentSolverDomain_)
                        )
                    );
                    localStartTime = tcTemp.localValue(0);
                }
                break;
            case mtsLatestTimeAllDomains:
                localStartTime = globalTime;
                break;
        }
    }

    startTime_ = localStartTime;

    globalTimeOffset_ = globalTime - startTime_;

    // Give multiDictRegistry a time value (required for regIOobject::write()
    // to case/[timeValue]
    multiDictRegistry_.setTime(startTime_, 0);

    // Copy the source data and any previous time directories to
    // case/[localTime]
    forAll(tcSource, i)
    {
        fileName copyMe(findInstancePath(tcSource, i));
        cp(copyMe, multiDictRegistry_.path());
    }
    mv
    (
        multiDictRegistry_.path()/sourcePath.name(),
        multiDictRegistry_.path()/Time::timeName(startTime_)
    );

    // If the source data was in a different domain, swap the boundary conditions
    if (sourcePath.path().path().name() != currentSolverDomain_)
    {
        swapBoundaryConditions
        (
            multiDictRegistry_.path()/Time::timeName(startTime_),
            currentSolverDomain_
        );
    }

    // Determine localEndTime and stopAtSetting
    word stopAtSetting(setLocalEndTime());

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
    // Start with timePrecision, it may affect writePrecision, which affects
    // all other settings
    if (multiSolverControl_.found("timePrecision"))
    {
        unsigned int newPrecision
        (
            readUint(multiSolverControl_.lookup("timePrecision"))
        );

        // Increase writePrecision if less than timePrecision, otherwise
        // resuming a run will fail
        if (IOstream::defaultPrecision() < (newPrecision + 1))
        {
            IOstream::defaultPrecision(newPrecision + 1);
        }

        newControlDict.set
        (
            "timePrecision",
            newPrecision
        );
    }

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

    // Write the dictionary to the case directory
    newControlDict.regIOobject::write();

    // Copy files in the superLoop/superLoopData directory of tcSource to
    // constant/superLoopData
    cp
    (
        sourcePath.path()/"superLoopData",
        multiDictRegistry_.path()/multiDictRegistry_.constant()
    );

    swapDictionaries(currentSolverDomain_);
    initialized_ = true;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
