    // Read field bounds
    dictionary fieldBounds = mesh.solutionDict().subDict("fieldBounds");

    // Pressure bounds
    dimensionedScalar pMin("pMin", p.dimensions(), SMALL);
    dimensionedScalar pMax("pMax", p.dimensions(), 1e10);
    fieldBounds.lookup(p.name()) >> pMin.value() >> pMax.value();

    // Temperature bounds
    dimensionedScalar TMin("TMin", T.dimensions(), SMALL);
    dimensionedScalar TMax("TMax", T.dimensions(), 1e10);
    fieldBounds.lookup(T.name()) >> TMin.value() >> TMax.value();

    // Density bounds
    dimensionedScalar rhoMin("rhoMin", rho.dimensions(), SMALL);
    dimensionedScalar rhoMax("rhoMax", rho.dimensions(), 1e10);
    fieldBounds.lookup(rho.name()) >> rhoMin.value() >> rhoMax.value();
