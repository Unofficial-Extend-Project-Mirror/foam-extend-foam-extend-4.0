    // Read field bounds
    dictionary fieldBounds = mesh.solutionDict().subDict("fieldBounds");

    // Pressure bounds
    dimensionedScalar pMin("pMin", p.dimensions(), 0);
    dimensionedScalar pMax("pMax", p.dimensions(), 0);

    fieldBounds.lookup(p.name()) >> pMin.value() >> pMax.value();

    // Velocity bound
    dimensionedScalar UMax("UMax", U.dimensions(), 0);

    fieldBounds.lookup(U.name()) >> UMax.value();
    dimensionedScalar smallU("smallU", dimVelocity, 1e-10);
