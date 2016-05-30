#   include "createTime.H"
#   include "createMesh.H"

    simpleControl simple(mesh);

#   include "createFields_scalarTransportFoam.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating scalar transport\n" << endl;

#   include "CourantNo.H"

    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        while (simple.correctNonOrthogonal())
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

