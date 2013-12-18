    // *** Adding data sources ***
        // Add time
        //  We could try addSource(runTime), but since runTime.name()
        //  constantly changes, we assign it our own name:
        eqns.scalarSources().addSource(runTime.value(), "t", dimTime);
        
        // Add mesh coordinates (cell centres) gives variable names:
        // C.x, C.y, C.z
        eqns.vectorSources().addSource(mesh.C());
        
        // Add mesh volumes
        eqns.scalarSources().addSource(mesh.V());
        
        // Add simpleFoam's existing variables
        eqns.scalarSources().addSource(p);
        eqns.vectorSources().addSource(U);
        
        // Adding a "derived" variable - one that exists only temporarily e.g.:
        // turbulence->R() - since it does not permanently exist, this won't
        // work:
        //      eqns.symmTensorSources().addSource(turbulence->R());
        // You have to create your own permanent variable, and update it at
        // every timestep:
        volSymmTensorField eqnSrcR
        (
            IOobject
            (
                "R_",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            turbulence->R()
        );
        eqns.symmTensorSources().addSource(eqnSrcR);

        // Add a dictionary source
        IOdictionary equationDict
        (
            IOobject
            (
                "equationDict",
                runTime.constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        eqns.addSource(equationDict);
        
        // Generic sources for demo purposes:
        
        // Scalars
        scalar sA(1.0);
        scalar sB(2.0);
        scalar sC(3.0);
        
        eqns.scalarSources().addSource(sA, "sA");
        eqns.scalarSources().addSource(sB, "sB");
        eqns.scalarSources().addSource(sC, "sC");
        
        // Dimensioned scalars
        dimensionedScalar dsA("dsA", dimless, 4.0);
        dimensionedScalar dsB("dsB", dimless, 5.0);
        dimensionedScalar dsC("dsC", dimless, 6.0);

        eqns.scalarSources().addSource(dsA);
        eqns.scalarSources().addSource(dsB);
        eqns.scalarSources().addSource(dsC);

        // scalarFields
        scalarField sfA(mesh.nCells(), 0.0);
        scalarField sfB(mesh.nCells(), 0.0);
        scalarField sfC(mesh.nCells(), 0.0);

        eqns.scalarSources().addSource(sfA, "sfA");
        eqns.scalarSources().addSource(sfB, "sfB");
        eqns.scalarSources().addSource(sfC, "sfC");

        // volScalarFields
        volScalarField vsfA
        (
            IOobject
            (
                "vsfA",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("vsfA", dimless, 0.0)
        );
        volScalarField vsfB
        (
            IOobject
            (
                "vsfB",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("vsfB", dimless, 0.0)
        );
        volScalarField vsfC
        (
            IOobject
            (
                "vsfC",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("vsfC", dimless, 0.0)
        );

        eqns.scalarSources().addSource(vsfA);
        eqns.scalarSources().addSource(vsfB);
        eqns.scalarSources().addSource(vsfC);

        // Vectors
        vector vA(100.0, 100.0, 100.0);
        vector vB(200.0, -200.0, 200.0);
        vector vC(300.0, 300.0, -300.0);
        eqns.vectorSources().addSource(vA, "vA");
        eqns.vectorSources().addSource(vB, "vB");
        eqns.vectorSources().addSource(vC, "vC");
        
        dimensionedVector dvA("dvA", dimless, vector(400.0, 400.0, 400.0));
        dimensionedVector dvB("dvB", dimless, vector(500.0, -500.0, 500.0));
        dimensionedVector dvC("dvC", dimless, vector(600.0, 600.0, -600.0));
        eqns.vectorSources().addSource(dvA);
        eqns.vectorSources().addSource(dvB);
        eqns.vectorSources().addSource(dvC);

        // vectorFields
        vectorField vfA(mesh.nCells(), vector(0.0, 0.0, 0.0));
        vectorField vfB(mesh.nCells(), vector(0.0, 0.0, 0.0));
        vectorField vfC(mesh.nCells(), vector(0.0, 0.0, 0.0));
        eqns.vectorSources().addSource(vfA, "vfA");
        eqns.vectorSources().addSource(vfB, "vfB");
        eqns.vectorSources().addSource(vfC, "vfC");

        volVectorField vvfA
        (
            IOobject
            (
                "vvfA",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedVector("vvfA", dimless, vector(0.0, 0.0, 0.0))
        );
        volVectorField vvfB
        (
            IOobject
            (
                "vvfB",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedVector("vvfB", dimless, vector(0.0, 0.0, 0.0))
        );
        volVectorField vvfC
        (
            IOobject
            (
                "vvfC",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedVector("vvfC", dimless, vector(0.0, 0.0, 0.0))
        );
        eqns.vectorSources().addSource(vvfA);
        eqns.vectorSources().addSource(vvfB);
        eqns.vectorSources().addSource(vvfC);

        // Tensor
        tensor tA
        (
            10000.0, 10000.0, 10000.0,
            10000.0, 10000.0, 10000.0,
            10000.0, 10000.0, 10000.0
        );
        tensor tB
        (
            20000.0, 20000.0, 20000.0,
            20000.0, 20000.0, 20000.0,
            20000.0, 20000.0, 20000.0
        );
        tensor tC
        (
            30000.0, 30000.0, 30000.0,
            30000.0, 30000.0, 30000.0,
            30000.0, 30000.0, 30000.0
        );
        eqns.tensorSources().addSource(tA, "tA");
        eqns.tensorSources().addSource(tB, "tB");
        eqns.tensorSources().addSource(tC, "tC");

        dimensionedTensor dtA
        (
            "dtA",
            dimless,
            tensor
            (
                40000.0, 40000.0, 40000.0,
                40000.0, 40000.0, 40000.0,
                40000.0, 40000.0, 40000.0
            )
        );
        dimensionedTensor dtB
        (
            "dtB",
            dimless,
            tensor
            (
                50000.0, 50000.0, 50000.0,
                50000.0, 50000.0, 50000.0,
                50000.0, 50000.0, 50000.0
            )
        );
        dimensionedTensor dtC
        (
            "dtC",
            dimless,
            tensor
            (
                60000.0, 60000.0, 60000.0,
                60000.0, 60000.0, 60000.0,
                60000.0, 60000.0, 60000.0
            )
        );
        eqns.tensorSources().addSource(dtA);
        eqns.tensorSources().addSource(dtB);
        eqns.tensorSources().addSource(dtC);

        // tensorFields
        tensorField tfA
        (
            mesh.nCells(),
            tensor
            (
                0.0, 0.0, 0.0,
                0.0, 0.0, 0.0,
                0.0, 0.0, 0.0
            )
        );
        tensorField tfB
        (
            mesh.nCells(),
            tensor
            (
                0.0, 0.0, 0.0,
                0.0, 0.0, 0.0,
                0.0, 0.0, 0.0
            )
        );
        tensorField tfC
        (
            mesh.nCells(),
            tensor
            (
                0.0, 0.0, 0.0,
                0.0, 0.0, 0.0,
                0.0, 0.0, 0.0
            )
        );
        eqns.tensorSources().addSource(tfA, "tfA");
        eqns.tensorSources().addSource(tfB, "tfB");
        eqns.tensorSources().addSource(tfC, "tfC");

        // volTensorFields
        volTensorField vtfA
        (
            IOobject
            (
                "vtfA",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedTensor
            (
                "vtfA",
                dimless,
                tensor
                (
                    0.0, 0.0, 0.0,
                    0.0, 0.0, 0.0,
                    0.0, 0.0, 0.0
                )
            )
        );
        volTensorField vtfB
        (
            IOobject
            (
                "vtfB",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedTensor
            (
                "vtfB",
                dimless,
                tensor
                (
                    0.0, 0.0, 0.0,
                    0.0, 0.0, 0.0,
                    0.0, 0.0, 0.0
                )
            )
        );
        volTensorField vtfC
        (
            IOobject
            (
                "vtfC",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedTensor
            (
                "vtfC",
                dimless,
                tensor
                (
                    0.0, 0.0, 0.0,
                    0.0, 0.0, 0.0,
                    0.0, 0.0, 0.0
                )
            )
        );
        eqns.tensorSources().addSource(vtfA);
        eqns.tensorSources().addSource(vtfB);
        eqns.tensorSources().addSource(vtfC);

        // And so on and so forth, including:
        //  - Types;
        //  - dimensionedTypes;
        //  - Fields;
        //  - DimensionedFields; and
        //  - GeometricFields,
        // where Types can be:
        //  - scalar
        //  - vector
        //  - tensor
        //  - diagTensor
        //  - symmTensor
        //  - sphericalTensor

    // *** Reading in the equations ***

        // You can read from any dictionary, but I'll just reuse the source
        // dictionary we created above.
        // scalar equations
        eqns.readEquation(equationDict, "sOut");
        eqns.readEquation(equationDict, "dsOut");
        eqns.readEquation(equationDict, "dsfOut");
        eqns.readEquation(equationDict, "volSfOut");
        
        // vector equations
        eqns.readEquation(equationDict, "vOut.x");
        eqns.readEquation(equationDict, "vOut.y");
        eqns.readEquation(equationDict, "vOut.z");
        eqns.readEquation(equationDict, "dvOut.x");
        eqns.readEquation(equationDict, "dvOut.y");
        eqns.readEquation(equationDict, "dvOut.z");
        eqns.readEquation(equationDict, "dvfOut.x");
        eqns.readEquation(equationDict, "dvfOut.y");
        eqns.readEquation(equationDict, "dvfOut.z");
        eqns.readEquation(equationDict, "volVfOut.x");
        eqns.readEquation(equationDict, "volVfOut.y");
        eqns.readEquation(equationDict, "volVfOut.z");
        
        // tensor equations
        eqns.readEquation(equationDict, "tOut.xx");
        eqns.readEquation(equationDict, "tOut.xy");
        eqns.readEquation(equationDict, "tOut.xz");
        eqns.readEquation(equationDict, "tOut.yx");
        eqns.readEquation(equationDict, "tOut.yy");
        eqns.readEquation(equationDict, "tOut.yz");
        eqns.readEquation(equationDict, "tOut.zx");
        eqns.readEquation(equationDict, "tOut.zy");
        eqns.readEquation(equationDict, "tOut.zz");
        eqns.readEquation(equationDict, "dtOut.xx");
        eqns.readEquation(equationDict, "dtOut.xy");
        eqns.readEquation(equationDict, "dtOut.xz");
        eqns.readEquation(equationDict, "dtOut.yx");
        eqns.readEquation(equationDict, "dtOut.yy");
        eqns.readEquation(equationDict, "dtOut.yz");
        eqns.readEquation(equationDict, "dtOut.zx");
        eqns.readEquation(equationDict, "dtOut.zy");
        eqns.readEquation(equationDict, "dtOut.zz");
        eqns.readEquation(equationDict, "dtfOut.xx");
        eqns.readEquation(equationDict, "dtfOut.xy");
        eqns.readEquation(equationDict, "dtfOut.xz");
        eqns.readEquation(equationDict, "dtfOut.yx");
        eqns.readEquation(equationDict, "dtfOut.yy");
        eqns.readEquation(equationDict, "dtfOut.yz");
        eqns.readEquation(equationDict, "dtfOut.zx");
        eqns.readEquation(equationDict, "dtfOut.zy");
        eqns.readEquation(equationDict, "dtfOut.zz");
        eqns.readEquation(equationDict, "volTfOut.xx");
        eqns.readEquation(equationDict, "volTfOut.xy");
        eqns.readEquation(equationDict, "volTfOut.xz");
        eqns.readEquation(equationDict, "volTfOut.yx");
        eqns.readEquation(equationDict, "volTfOut.yy");
        eqns.readEquation(equationDict, "volTfOut.yz");
        eqns.readEquation(equationDict, "volTfOut.zx");
        eqns.readEquation(equationDict, "volTfOut.zy");
        eqns.readEquation(equationDict, "volTfOut.zz");

    // *** Create output objects ***

        // output dictionary (for objects that don't output themselves)
        IOdictionary eqnOutputDict
        (
            IOobject
            (
                "eqnOutput",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            )
        );
    
        // scalars
        scalar scalarOut;
        dimensionedScalar dimensionedScalarOut
        (
            "dsOut",
            eqns.evaluateDimensions("dsOut"),
            0.0
        );

        DimensionedField<scalar, volMesh> dimensionedScalarFieldOut
        (
            IOobject
            (
                "dsfOut",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar
            (
                "dsfOut",
                eqns.evaluateDimensions("dsfOut"),
                0.0
            )
        );

        volScalarField volScalarFieldOut
        (
            IOobject
            (
                "volSfOut",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar
            (
                "volSfOut",
                eqns.evaluateDimensions("volSfOut"),
                0.0
            )
        );

        // vector
        vector vectorOut;
        dimensionedVector dimensionedVectorOut
        (
            "dvOut",
            eqns.evaluateDimensions("dvOut.x"),
            vector(0.0, 0.0, 0.0)
        );

        DimensionedField<vector, volMesh> dimensionedVectorFieldOut
        (
            IOobject
            (
                "dvfOut",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedVector
            (
                "dvfOut",
                eqns.evaluateDimensions("dvfOut.x"),
                vector(0.0, 0.0, 0.0)
            )
        );

        volVectorField volVectorFieldOut
        (
            IOobject
            (
                "volVfOut",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedVector
            (
                "volVfOut",
                eqns.evaluateDimensions("volVfOut.x"),
                vector(0.0, 0.0, 0.0)
            )
        );

        // tensors
        tensor tensorOut;
        dimensionedTensor dimensionedTensorOut
        (
            "dtOut",
            eqns.evaluateDimensions("dtOut.xx"),
            tensor
            (
                0.0, 0.0, 0.0,
                0.0, 0.0, 0.0,
                0.0, 0.0, 0.0
            )
        );

        DimensionedField<tensor, volMesh> dimensionedTensorFieldOut
        (
            IOobject
            (
                "dtfOut",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedTensor
            (
                "dtfOut",
                eqns.evaluateDimensions("dtfOut.xx"),
                tensor
                (
                    0.0, 0.0, 0.0,
                    0.0, 0.0, 0.0,
                    0.0, 0.0, 0.0
                )
            )
        );
        volTensorField volTensorFieldOut
        (
            IOobject
            (
                "volTfOut",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedTensor
            (
                "volTfOut",
                eqns.evaluateDimensions("volTfOut.xx"),
                tensor
                (
                    0.0, 0.0, 0.0,
                    0.0, 0.0, 0.0,
                    0.0, 0.0, 0.0
                )
            )
        );

