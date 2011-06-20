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


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::multiSolver::setMultiSolverControls()
{
    initialStartFrom_ = misLatestTime;
    if (multiSolverControl_.found("initialStartFrom"))
    {
        initialStartFrom_ = initialStartFromControlsNames_.read
        (
            multiSolverControl_.lookup("initialStartFrom")
        );
    }
    
    if (multiSolverControl_.found("startTime"))
    {
        initialStartTime_ = readScalar(multiSolverControl_.lookup("startTime"));
        if (initialStartTime_ < 0)
        {
            FatalErrorIn("multiSolver::setMultiSolverControls")
                << "'startTime' in multiControlDict/multiSolverControl cannot "
                << "be negative."
                << abort(FatalError);
        }
    }
    else if
    (
        (initialStartFrom_ == misStartTime)
     || (initialStartFrom_ == misStartTimeInStartDomain)
     || (initialStartFrom_ == misStartTimeInStartDomainInStartSuperLoop)
    )
    {
        FatalIOErrorIn
        (
            "multiSolver::setMultiSolverControls", multiSolverControl_
        )
            << "'startTime' is required in multiControlDict/multiSolverControl "
            << "if 'initialStartFrom' is set to 'startTime', "
            << "'startTimeInStartDomain', or "
            << "'startTimeInStartDomainInStartStuperLoop'"
            << exit(FatalIOError);
    }

    if (multiSolverControl_.found("startDomain"))
    {
        startDomain_ = word(multiSolverControl_.lookup("startDomain"));
    }
    else if 
    (
        (initialStartFrom_ == misFirstTimeInStartDomain)
     || (initialStartFrom_ == misFirstTimeInStartDomainInStartSuperLoop)
     || (initialStartFrom_ == misStartTimeInStartDomain)
     || (initialStartFrom_ == misStartTimeInStartDomainInStartSuperLoop)
     || (initialStartFrom_ == misLatestTimeInStartDomain)
     || (initialStartFrom_ == misLatestTimeInStartDomainInStartSuperLoop)
    )
    {
        FatalIOErrorIn("multiSolver::setMultiSolverControls", multiSolverControl_)
            << "'startDomain' is required in "
            << "multiControlDict/multiSolverControl if 'initialStartFrom' is "
            << "set to 'firstTimeInStartDomain', "
            << "'firstTimeInStartDomainInStartSuperLoop', "
            << "'startTimeInStartDomain', "
            << "'startTimeInStartDomainInStartSuperLoop', "
            << "'latestTimeInStartDomain', or "
            << "'latestTimeInStartDomainInStartSuperLoop'."
            << abort(FatalError);
    }
    
    finalStopAt_ = mfsEndTime;
    if (multiSolverControl_.found("finalStopAt"))
    {
        finalStopAt_ = finalStopAtControlsNames_.read
        (
            multiSolverControl_.lookup("finalStopAt")
        );
    }
    
    if (multiSolverControl_.found("endDomain"))
    {
        endDomain_ = word(multiSolverControl_.lookup("endDomain"));
    }
    else if
    (
        (finalStopAt_ == mfsEndTimeInEndDomain)
     || (finalStopAt_ == mfsEndTimeInEndDomainInEndSuperLoop)
    )
    {
        FatalErrorIn("multiSolver::setMultiSolverControls")
            << "endTime is required in multiControlDict/multiSolverControl if "
            << "finalStopAt is set to 'endTimeInEndDomain', or "
            << "'endTimeInEndDomainInEndSuperLoop'."
            << abort(FatalError);
    }

    if (multiSolverControl_.found("endTime"))
    {
        finalEndTime_ =
            readScalar(multiSolverControl_.lookup("endTime"));
    }
    else if
    (
        (finalStopAt_ == mfsEndTime)
     || (finalStopAt_ == mfsEndTimeInEndDomain)
     || (finalStopAt_ == mfsEndTimeInEndDomainInEndSuperLoop)
    )
    {
        FatalErrorIn("multiSolver::setMultiSolverControls")
            << "'endTime' is required in "
            << "multiControlDict/multiSolverControl if 'finalStopAt' is set to "
            << "'endTime', 'endTimeInEndDomain', or "
            << "'endTimeInEndDomainInEndSuperLoop'."
            << abort(FatalError);
    }

    if (multiSolverControl_.found("startSuperLoop"))
    {
        startSuperLoop_ =
            readLabel(multiSolverControl_.lookup("startSuperLoop"));
    }
    else if
    (
        (initialStartFrom_ == misFirstTimeInStartDomainInStartSuperLoop)
     || (initialStartFrom_ == misStartTimeInStartDomainInStartSuperLoop)
     || (initialStartFrom_ == misLatestTimeInStartDomainInStartSuperLoop)
    )
    {
        FatalErrorIn("multiSolver::setMultiSolverControls")
            << "'startSuperLoop' is required in "
            << "multiControlDict/multiSolverControl if 'initialStartFrom' is "
            << "set to 'firstTimeInStartDomainInSuperLoop', "
            << "'startTimeInStartDomainInStartSuperLoop', or "
            << "'latestTimeInStartDomainInStartSuperLoop'."
            << abort(FatalError);
    }
    
    if (multiSolverControl_.found("endSuperLoop"))
    {
        endSuperLoop_ = 
            readLabel(multiSolverControl_.lookup("endSuperLoop"));
    }
    else if
    (
        (finalStopAt_ == mfsEndSuperLoop)
     || (finalStopAt_ == mfsEndTimeInEndDomainInEndSuperLoop)
    )
    {
        FatalErrorIn("multiSolver::setMultiSolverControls")
            << "'endSuperLoops' is required in "
            << "multiControlDict/multiSolverControl if 'finalStopAt' is set to "
            << "'endSuperLoop' or 'endTimeInEndDomainInEndSuperLoop'."
            << abort(FatalError);
    }

    multiDictsRunTimeModifiable_ = true;
    if (multiSolverControl_.found("multiDictsRunTimeModifiable"))
    {
        multiDictsRunTimeModifiable_ =
            multiSolverControl_.lookup("multiDictsRunTimeModifiable");
    }

    prefixes_.clear();
    prefixes_ = solverDomains_.toc();
    if
    (
        (prefixes_.size() == 0)
     || ((prefixes_.size() == 1) && (prefixes_[0] == "default"))
    )
    {
        FatalErrorIn("multiSolver::setMultiSolverControls")
            << "No solver domains found in multiControlDict.  Expecting "
            << "subdictionary solverDomains to contain at least one entry "
            << "other than 'default'."
            << abort(FatalError);
    }
    
    dictionary solverDomainsDefault
    (
        solverDomains_.found("default")
      ? solverDomains_.subDict("default")
      : dictionary()
    );
}


void Foam::multiSolver::setSolverDomainControls(const word& solverDomainName)
{
    currentSolverDomainDict_.clear();
    if (solverDomains_.found("default"))
    {
        currentSolverDomainDict_.merge(solverDomains_.subDict("default"));
    }
    currentSolverDomainDict_.merge(solverDomains_.subDict(solverDomainName));

    startFrom_ = mtsLatestTimeAllDomains;
    if (currentSolverDomainDict_.found("startFrom"))
    {
        startFrom_ = startFromControlsNames_.read
        (
            currentSolverDomainDict_.lookup("startFrom")
        );
    }

    if (currentSolverDomainDict_.found("startTime"))
    {
        startTime_ = readScalar(currentSolverDomainDict_.lookup("startTime"));
        if (startTime_ < 0)
        {
            FatalErrorIn("multiSolver::setSolverDomainControls")
                << "'startTime' in multiControlDict/solverDomains/" 
                << solverDomainName << " cannot be negative."
                << abort(FatalError);
        }
    }
    else if (startFrom_ == mtsStartTime)
    {
        FatalErrorIn("multiSolver::setSolverDomainControls")
            << "'startTime' not defined in solverDomain '" << solverDomainName
            << "' or 'default'.  startTime is required when startFrom "
            << "is set to 'startTime'."
            << abort(FatalError);
    }
    
    stopAt_ = msaEndTime;
    if (currentSolverDomainDict_.found("stopAt"))
    {
        stopAt_ = stopAtControlsNames_.read
        (
            currentSolverDomainDict_.lookup("stopAt")
        );
        if
        (
            (stopAt_ == msaIterations)
         && (currentSolverDomainDict_.found("adjustTimeStep"))
         && (
                readBool(currentSolverDomainDict_.lookup
                ("adjustTimeStep")) == true
            )
        )
        {
        FatalErrorIn("multiSolver::setSolverDomainControls")
            << "'stopAt' in multiControlDict/sovlerDomains cannot be set "
            << "to 'iterations' when 'adjustTimeStep' is 'true'."
            << abort(FatalError);
        }
    }
    
    endTime_ = 0;
    if (currentSolverDomainDict_.found("endTime"))
    {
        endTime_ = readScalar(currentSolverDomainDict_.lookup("endTime"));
    }

    if (currentSolverDomainDict_.found("iterations"))
    {
        iterations_ = readLabel
        (
            currentSolverDomainDict_.lookup("iterations")
        );
    }
    else if (stopAt_ == msaIterations)
    {
        FatalErrorIn("multiSolver::setSolverDomainControls")
            << "'iterations' not defined in solverDomain '"
            << solverDomainName << "' or 'default'.  iterations is required "
            << "when stopAt is set to iterations."
            << abort(FatalError);
    }
    
    if (currentSolverDomainDict_.found("elapsedTime"))
    {
        elapsedTime_ = readScalar(currentSolverDomainDict_.lookup("elapsedTime"));
    }
    else if (stopAt_ == msaElapsedTime)
    {
        FatalErrorIn("multiSolver::setSolverDomainControls")
            << "'elapsedTime' not defined in solverDomain '"
            << solverDomainName << "' or 'default'.  elapsedTime is required "
            << "when stopAt is set to elapsedTime."
            << abort(FatalError);
    }

    if (currentSolverDomainDict_.found("storeFields"))
    {
        storeFields_ = wordList
        (
            currentSolverDomainDict_.lookup("storeFields")
        );
    }
    
    purgeWriteSuperLoops_ = 0;
    if (currentSolverDomainDict_.found("purgeWriteSuperLoops"))
    {
        purgeWriteSuperLoops_ = readLabel
        (
            currentSolverDomainDict_.lookup("purgeWriteSuperLoops")
        );
    }
    
    if (currentSolverDomainDict_.found("deltaT"))
    {
        deltaT_ = readScalar
        (
            currentSolverDomainDict_.lookup("deltaT")
        );
    }
    else
    {
        FatalErrorIn("multiSolver::setSolverDomainControls")
            << "'deltaT' not defined in solverDomain '"
            << solverDomainName << "' or 'default'.  deltaT is required."
            << abort(FatalError);
    }
    
    if
    (
        currentSolverDomainDict_.found("timeFormat")
     || currentSolverDomainDict_.found("timePrecision")
    )
    {
        WarningIn("multiSolver::setSolverDomainControls")
            << "Dictionary entry 'timeFormat' or 'timePrecision' found in "
            << "multiControlDict/solverDomain subdictionaries and will be "
            << "ignored. This setting must be applied universally in "
            << "multiControlDict/multiSolverControl."
            << endl;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

