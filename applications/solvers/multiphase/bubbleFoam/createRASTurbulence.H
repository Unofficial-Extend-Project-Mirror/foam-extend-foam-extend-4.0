    IOdictionary RASProperties
    (
        IOobject
        (
            "RASProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );


    Switch turbulence
    (
        RASProperties.lookup("turbulence")
    );

    dictionary kEpsilonCoeffs
    (
        RASProperties.subDict("kEpsilonCoeffs")
    );


    scalar Cmu
    (
        readScalar(kEpsilonCoeffs.lookup("Cmu"))
    );

    scalar C1
    (
        readScalar(kEpsilonCoeffs.lookup("C1"))
    );

    scalar C2
    (
        readScalar(kEpsilonCoeffs.lookup("C2"))
    );

    scalar alphak
    (
        readScalar(kEpsilonCoeffs.lookup("alphak"))
    );

    scalar alphaEps
    (
        readScalar(kEpsilonCoeffs.lookup("alphaEps"))
    );

    dictionary wallFunctionCoeffs
    (
        RASProperties.subDict("wallFunctionCoeffs")
    );

    scalar kappa
    (
        readScalar(wallFunctionCoeffs.lookup("kappa"))
    );

    scalar E
    (
        readScalar(wallFunctionCoeffs.lookup("E"))
    );

    nearWallDist y(mesh);


    Info<< "Reading field k\n" << endl;
    volScalarField k
    (
        IOobject
        (
            "k",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field epsilon\n" << endl;
    volScalarField epsilon
    (
        IOobject
        (
            "epsilon",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    Info<< "Calculating field nutb\n" << endl;
    volScalarField nutb
    (
        IOobject
        (
            "nutb",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        Cmu*sqr(k)/epsilon
    );

    Info<< "Calculating field nuEffa\n" << endl;
    volScalarField nuEffa
    (
        IOobject
        (
            "nuEffa",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        sqr(Ct)*nutb + nua
    );

    Info<< "Calculating field nuEffb\n" << endl;
    volScalarField nuEffb
    (
        IOobject
        (
            "nuEffb",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        nutb + nub
    );
