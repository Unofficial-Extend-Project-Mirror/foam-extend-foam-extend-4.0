{
    rho = thermo.rho();

    // Thermodynamic density needs to be updated by psi*d(p) after the
    // pressure solution - done in 2 parts. Part 1:
    thermo.rho() -= psi*p;

    volScalarField rAU = 1.0/UEqn.A();
    U = rAU*UEqn.H();

    if (pZones.size() > 0)
    {
        // ddtPhiCorr not well defined for cases with porosity
        phi = fvc::interpolate(rho)*(fvc::interpolate(U) & mesh.Sf());
    }
    else
    {
        phi =
            fvc::interpolate(rho)
           *(
                (fvc::interpolate(U) & mesh.Sf())
              + fvc::ddtPhiCorr(rAU, rho, U, phi)
            );
    }

    while (piso.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvc::ddt(rho) + psi*correction(fvm::ddt(p))
          + fvc::div(phi)
          - fvm::laplacian(rho*rAU, p)
        ==
            parcels.Srho()
          + massSource.SuTot()
        );

        pEqn.solve(mesh.solutionDict().solver(p.select(piso.finalInnerIter())));

        if (piso.finalNonOrthogonalIter())
        {
            phi += pEqn.flux();
        }
    }

    // Second part of thermodynamic density update
    thermo.rho() += psi*p;

    #include "rhoEqn.H"
    #include "compressibleContinuityErrs.H"

    U -= rAU*fvc::grad(p);
    U.correctBoundaryConditions();
}
