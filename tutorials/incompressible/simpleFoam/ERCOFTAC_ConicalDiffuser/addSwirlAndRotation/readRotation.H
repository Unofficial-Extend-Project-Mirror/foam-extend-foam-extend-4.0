 
    Info << "\nReading swirlAndRotationProperties dictionary\n" << endl;
    IOdictionary swirlAndRotationProperties
    (
         IOobject
         (
            "swirlAndRotationProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
          )
    );

    Switch modifyBoundaries
    (
        swirlAndRotationProperties.lookup("modifyBoundaries")
    );

    Switch modifyInterior
    (
        swirlAndRotationProperties.lookup("modifyInterior")
    );

    dimensionedVector omega
    (
        swirlAndRotationProperties.lookup("omega")
    );

    dimensionedVector center
    (
        swirlAndRotationProperties.lookup("center")
    );

    wordList rotPatches(swirlAndRotationProperties.lookup("rotatingPatches"));

