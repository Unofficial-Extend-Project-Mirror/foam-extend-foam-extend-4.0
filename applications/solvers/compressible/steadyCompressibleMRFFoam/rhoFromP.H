{
    // Calculate density from pressure
    rho.storePrevIter();
    rho = thermo.rho()();

    // Bound rho
    volScalarField R = thermo.Cp() - thermo.Cv();

    volScalarField rhoMin = pMin/(R*TMax);
    volScalarField rhoMax = pMax/(R*TMin);

    rho = Foam::min(rho, rhoMax);
    rho = Foam::max(rho, rhoMin);
    rho.relax();
    rho.correctBoundaryConditions();
}
