{
    // Stabilisation of friction factor calculation
    // Friction factor is defined with standard gravity
    frictionFactor.internalField() =
        mag(2*9.81*sqr(manningField.internalField())/
        pow(mag(h.internalField()) + 1e-7, 1.0/3.0));
}
