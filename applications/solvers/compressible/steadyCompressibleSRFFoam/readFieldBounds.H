    // Read field bounds
    dictionary fieldBounds = mesh.solutionDict().subDict("fieldBounds");

    // Pressure bounds
    dimensionedScalar pMin("pMin", dimPressure, 0);
    dimensionedScalar pMax("pMax", dimPressure, GREAT);

    fieldBounds.lookup("p") >> pMin.value() >> pMax.value();

    // Temperature bounds
    dimensionedScalar TMin("TMin", dimTemperature, 0);
    dimensionedScalar TMax("TMax", dimTemperature, GREAT);

    fieldBounds.lookup("T") >> TMin.value() >> TMax.value();

    // Velocity bound
    dimensionedScalar UrelMax("UrelMax", dimVelocity, GREAT);

    fieldBounds.lookup(Urel.name()) >> UrelMax.value();
    dimensionedScalar smallUrel("smallUrel", dimVelocity, 1e-10);
