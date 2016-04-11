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
    dimensionedScalar UMax("UMax", U.dimensions(), GREAT);

    fieldBounds.lookup(U.name()) >> UMax.value();
    dimensionedScalar smallU("smallU", dimVelocity, 1e-10);
