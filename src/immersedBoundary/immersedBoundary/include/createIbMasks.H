//
// createIbMask.H
// ~~~~~~~~~~~~

    Info<< "Create immersed boundary cell mask" << endl;

    volScalarField cellIbMask
    (
        IOobject
        (
            "cellIbMask",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("one", dimless, 1)
    );

    Info<< "Create immersed boundary face mask" << endl;
    surfaceScalarField faceIbMask
    (
        IOobject
        (
            "faceIbMask",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("one", dimless, 1)
    );

    // Multiple IB patches currently not supported.
    // This required IB geometry update to re-combine the corrected fields in
    // Geometry correction functions in immersedBoundaryFvPatch
    // HJ, 30/Nov/2017
    {
        label nIbPatches = 0;

        forAll (mesh.boundary(), patchI)
        {
            if (isA<immersedBoundaryFvPatch>(mesh.boundary()[patchI]))
            {
                Info<< "Found immersed boundary patch " << patchI
                    << " named " << mesh.boundary()[patchI].name()
                    << endl;

                nIbPatches++;
                
                const immersedBoundaryFvPatch& ibPatch =
                    refCast<const immersedBoundaryFvPatch>
                    (
                        mesh.boundary()[patchI]
                    );

                cellIbMask = Foam::min(cellIbMask, ibPatch.gamma());
                faceIbMask *= Foam::min(faceIbMask, ibPatch.sGamma());
            }
        }

        if (nIbPatches > 1)
        {
            FatalErrorIn(args.executable())
                << "Multiple Immersed boundary patches currently not supported"
                << abort(FatalError);
        }
    }

    // Evaluate boundary conditions for IB masks
    cellIbMask.boundaryField().evaluateCoupled();
