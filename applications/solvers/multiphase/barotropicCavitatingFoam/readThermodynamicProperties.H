    Info<< "Reading thermodynamicProperties\n" << endl;

    IOdictionary thermodynamicProperties
    (
        IOobject
        (
            "thermodynamicProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar rhogSat
    (
        thermodynamicProperties.lookup("rhogSat")
    );

    dimensionedScalar rholSat
    (
        thermodynamicProperties.lookup("rholSat")
    );

    dimensionedScalar pSat
    (
        thermodynamicProperties.lookup("pSat")
    );

    dimensionedScalar psig
    (
        thermodynamicProperties.lookup("psig")
    );

    dimensionedScalar psil
    (
        thermodynamicProperties.lookup("psil")
    );
