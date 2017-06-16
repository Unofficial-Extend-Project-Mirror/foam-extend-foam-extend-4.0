    label pRefCell = 0;
    scalar pRefValue = 0;

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
