    label pRefCell = 0;
    scalar pRefValue = 0;

    // Number of outer correctors
    const label nOuterCorrectors = readLabel
    (
        mesh.solutionDict().subDict("blockSolver").lookup("nOuterCorrectors")
    );

    setRefCell
    (
        p,
        mesh.solutionDict().subDict("blockSolver"),
        pRefCell,
        pRefValue
    );

    mesh.solutionDict().subDict("blockSolver").readIfPresent
    (
        "convergence",
        convergenceCriterion
    );

    mesh.solutionDict().subDict("blockSolver").readIfPresent
    (
        "convergence",
        convergenceCriterion
    );

    mesh.solutionDict().subDict("blockSolver").readIfPresent
    (
        "convergence",
        convergenceCriterion
    );
