{
    // Decoupled patches
#   include "attachPatches.H"

    // Solid side
#   include "readSolidControls.H"

    for (int nonOrth = 0; nonOrth <= nNonOrthCorr; nonOrth++)
    {
        coupledFvScalarMatrix TEqns(2);

        // Add fluid equation
        TEqns.set
        (
            0,
            new fvScalarMatrix
            (
                fvm::ddt(T)
              + fvm::div(phi, T)
              - fvm::laplacian(DT, T)
            )
        );

        // Add solid equation
        TEqns.set
        (
            1,
            new fvScalarMatrix
            (
                fvm::ddt(Tsolid) - fvm::laplacian(DTsolid, Tsolid)
            )
        );

        TEqns.solve();
    }
}
