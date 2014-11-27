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
    dimensionedScalar UMax("UMax", dimVelocity, GREAT);

    fieldBounds.lookup(U.name()) >> UMax.value();
    dimensionedScalar smallU("smallU", dimVelocity, 1e-10);
