// Reading of the controling parameters of the Picard loop and of the
// Stabilisation cycles loop for the adaptive time step procedure.
// Localisation of the data within the considered case: system/fvSolution
const dictionary& picardDict = mesh.solutionDict().subDict("Picard");

// Maximum number of Picard iterations.
const int nIterPicard =
    picardDict.lookupOrDefault<int>("nIterPicard", 1);

// Maximum number of stabilisation cycles.
const int nMaxCycle =
    picardDict.lookupOrDefault<int>("nMaxCycle", 1);

// Number of time iterations with low number of Picard iteration beyond which
// the time step is increased.
const int stabilisationThreshold =
    picardDict.lookupOrDefault<int>("stabilisationThreshold", 1);

// Exit criterion for the Picard loop.
const double precPicard =
    picardDict.lookupOrDefault<double>("precPicard", 1.);

// Factor of increase/decrease of the time step.
const double tFact =
    picardDict.lookupOrDefault<double>("tFact", 1.);
