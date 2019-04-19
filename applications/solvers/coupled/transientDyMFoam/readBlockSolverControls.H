{
    dictionary blockSolverDict = mesh.solutionDict().subDict("blockSolver");

    setRefCell
    (
        p,
        blockSolverDict,
        pRefCell,
        pRefValue
    );

    blockSolverDict.readIfPresent
    (
        "nOuterCorrectors",
        nOuterCorrectors
    );

    blockSolverDict.readIfPresent
    (
        "correctPhi",
        correctPhi
    );

    blockSolverDict.readIfPresent
    (
        "checkMeshCourantNo",
        checkMeshCourantNo
    );
}
