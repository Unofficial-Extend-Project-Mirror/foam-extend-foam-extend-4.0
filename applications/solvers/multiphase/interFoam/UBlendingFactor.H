    surfaceScalarField gammaf = fvc::interpolate(gamma);
    surfaceScalarField UBlendingFactor
    (
        "UBlendingFactor",
        sqrt(max(min(4*gammaf*(1.0 - gammaf), 1.0), 0.0))
    );
