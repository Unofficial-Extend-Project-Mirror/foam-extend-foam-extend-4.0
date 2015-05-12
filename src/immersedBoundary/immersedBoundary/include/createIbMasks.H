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

    volScalarField cellIbMaskExt
    (
        IOobject
        (
            "cellIbMaskExt",
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

    forAll (mesh.boundary(), patchI)
    {
        if (isA<immersedBoundaryFvPatch>(mesh.boundary()[patchI]))
        {
            Info<< "Found immersed boundary patch " << patchI
                << " named " << mesh.boundary()[patchI].name()
                << endl;

            const immersedBoundaryFvPatch& ibPatch =
                refCast<const immersedBoundaryFvPatch>
                (
                    mesh.boundary()[patchI]
                );

            cellIbMask *= ibPatch.gamma();
            cellIbMaskExt *= ibPatch.gammaExt();
            faceIbMask *= ibPatch.sGamma();
        }
    }

    // Evaluate boundary conditions for IB masks
    cellIbMask.boundaryField().evaluateCoupled();
    cellIbMaskExt.boundaryField().evaluateCoupled();
