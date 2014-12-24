//
// createEquationReader.H
// ~~~~~~~~~~~~~~~~~~~~~~

    Foam::IOEquationReader eqns
    (
        IOobject
        (
            "equations",
            runTime.timeName(),
            runTime,
            IOobject::NO_READ,
            runTime.controlDict().found("outputEquations")
          ? Switch(runTime.controlDict().lookup("outputEquations")) == true
            ? IOobject::AUTO_WRITE
            : IOobject::NO_WRITE
          : IOobject::NO_WRITE
        ),
        runTime.controlDict().found("outputEquationDataSources")
      ? bool(Switch(runTime.controlDict().lookup("outputEquationDataSources")))
      : false
    );
