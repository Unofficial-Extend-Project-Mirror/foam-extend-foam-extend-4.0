    // Read field bounds
    dictionary fieldBounds = mesh.solutionDict().subDict("fieldBounds");

    // Pressure bounds
    dimensionedScalar pMin("pMin", p.dimensions(), 0);
    dimensionedScalar pMax("pMax", p.dimensions(), GREAT);

    fieldBounds.lookup(p.name()) >> pMin.value() >> pMax.value();

    // Temperature bounds
    dimensionedScalar TMin("TMin", T.dimensions(), 0);
    dimensionedScalar TMax("TMax", T.dimensions(), GREAT);

    fieldBounds.lookup(T.name()) >> TMin.value() >> TMax.value();

    // Velocity bound
    dimensionedScalar UrelMax("UrelMax", Urel.dimensions(), GREAT);

    fieldBounds.lookup(Urel.name()) >> UrelMax.value();
    dimensionedScalar smallUrel("smallUrel", dimVelocity, 1e-10);
