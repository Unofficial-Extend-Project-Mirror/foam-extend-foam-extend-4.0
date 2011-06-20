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
    timeCluster tcSource;
    switch (initialStartFrom_)
    {
        case misFirstTime:
            tcSource = findClosestGlobalTime
            (
                0, readSuperLoopTimes(currentSolverDomain_, -1)
            );
            break;
        case misFirstTimeInStartDomain:
            tcSource = findClosestGlobalTime
            (
                0, readSuperLoopTimes(startDomain_, -1)
            );
            break;
        case misFirstTimeInStartDomainInStartSuperLoop:
            tcSource = findClosestGlobalTime
            (
                0, readSuperLoopTimes(startDomain_, startSuperLoop_)
            );
            break;
        case misStartTime:
            if (initialStartTime_ == 0)
            {
                tcSource = findClosestGlobalTime
                (
                    initialStartTime_,
                    readSuperLoopTimes(currentSolverDomain_, -1)
                );
            }
            else
            {
                tcSource = findClosestGlobalTime
                (
                    initialStartTime_, readAllTimes()
                );
            }
            break;
        case misStartTimeInStartDomain:
            tcSource = findClosestLocalTime
            (
                initialStartTime_, readSolverDomainTimes(startDomain_)
            );
            break;
        case misStartTimeInStartDomainInStartSuperLoop:
            tcSource = findClosestLocalTime
            (
                initialStartTime_,
                readSuperLoopTimes(startDomain_, startSuperLoop_)
            );
            break;
        case misLatestTime:
            tcSource = findLatestGlobalTime(readAllTimes());
            break;
        case misLatestTimeInStartDomain:
            tcSource = findLatestLocalTime(readSolverDomainTimes(startDomain_));
            break;
        case misLatestTimeInStartDomainInStartSuperLoop:
            tcSource = findLatestLocalTime
            (
                readSuperLoopTimes(startDomain_, startSuperLoop_)
            );
            break;
    }

    if (!tcSource.times().size())
    {
        // No relevant data found, set to initial conditions
        tcSource = timeCluster
        (
            Time::findTimes
            (
                multiDictRegistry_.path()/"multiSolver"/currentSolverDomain_
                    /"initial"
            ),
            0,
            -1, // superLoop of -1 signifies "initial" directory
            currentSolverDomain_
        );
    }

    fileName sourcePath(findInstancePath(tcSource, 0));
    superLoop_ = tcSource.superLoop();
    // If starting from initial conditions, superLoop_ = -1
    if (superLoop_ < 0) superLoop_ = 0;
    scalar globalTime(tcSource.globalValue(0));
    scalar localStartTime(tcSource.localValue(0));

    // Now to apply the exceptions if currentSolverDomain_ != data source
    // solverDomain (see long comment above).
    if (sourcePath.path().path().name() != currentSolverDomain_)
    {
        superLoop_++;
        
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

    // Copy the source data to case/[localTime]
    cp(sourcePath, multiDictRegistry_.path());
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

