//
// createIbMask.H
// ~~~~~~~~~~~~

    // Multiple IB patches currently not supported.
    // This required IB geometry update to re-combine the corrected fields in
    // Geometry correction functions in immersedBoundaryFvPatch
    // HJ, 30/Nov/2017

    Info<< "Create immersed boundary cell mask" << endl;

    volScalarField gamma
    (
        IOobject
        (
            "gamma",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("one", dimless, 1)
    );

    Info<< "Create immersed boundary face mask" << endl;
    surfaceScalarField sGamma
    (
        IOobject
        (
            "sGamma",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("one", dimless, 1)
    );

#   include "updateIbMasks.H"
