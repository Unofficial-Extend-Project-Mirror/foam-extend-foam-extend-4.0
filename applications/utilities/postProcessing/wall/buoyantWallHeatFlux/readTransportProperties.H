    singlePhaseTransportModel laminarTransport(U, phi);

    // Laminar Prandtl number
    const dimensionedScalar Pr(laminarTransport.lookup("Pr"));

    // Turbulent Prandtl number
    const dimensionedScalar Prt(laminarTransport.lookup("Prt"));

    // Reference density
    const dimensionedScalar rhoRef(laminarTransport.lookup("rhoRef"));

    // Specific heat capacity
    const dimensionedScalar c(laminarTransport.lookup("c"));
