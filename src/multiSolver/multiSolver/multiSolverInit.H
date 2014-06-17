// Common elements in the constructor initialization list for multiSolver

    multiDicts_(0),

    multiSolverControl_(multiControlDict_.subDict("multiSolverControl")),
    solverDomains_(multiControlDict_.subDict("solverDomains")),
    currentSolverDomainDict_(),
    prefixes_(solverDomains_.toc()),
    initialized_(false),
    forcedEnd_(false),
    noSaveSinceSuperLoopIncrement_(false),
    manageLocalRoot_(false),

    superLoop_(0),
    currentSolverDomain_("default"),

    initialStartFrom_(misLatestTime),
    initialStartTime_(0),
    startDomain_("default"),
    startSuperLoop_(0),
    finalStopAt_(mfsEndTime),
//    endDomain_,
    finalEndTime_(0),

    endSuperLoop_(0),

    multiDictsRunTimeModifiable_(true),

    globalTimeOffset_(0),
    globalIndex_(0),

    purgeWriteSuperLoops_(prefixes_.size()),

//    timeValueStartFrom_,
//    scalar startTime_,
//    stopAtControls stopAt_,
//    scalar endTime_,
//    label endIterations_

    storeFields_(0)

