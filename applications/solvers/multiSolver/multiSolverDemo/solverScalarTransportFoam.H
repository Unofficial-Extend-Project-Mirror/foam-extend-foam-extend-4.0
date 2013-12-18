#   include "createTime.H"
#   include "createMesh.H"
#   include "createFields_scalarTransportFoam.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating scalar transport\n" << endl;

#   include "CourantNo.H"

    for (runTime++; !runTime.end(); runTime++)
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

    dictionary simple = mesh.solutionDict().subDict("SIMPLE");

    int nNonOrthCorr =
        simple.lookupOrDefault<int>("nNonOrthogonalCorrectors", 0);

    bool momentumPredictor =
        simple.lookupOrDefault<Switch>("momentumPredictor", true);

        for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
        {
            solve
            (
                fvm::ddt(T)
              + fvm::div(phi, T)
              - fvm::laplacian(DT, T)
            );
        }

        runTime.write();
    }

    Info<< "End\n" << endl;

