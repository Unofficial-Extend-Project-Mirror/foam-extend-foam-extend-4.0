{
    volScalarField rUA = 1.0/UEqn.A();
    surfaceScalarField rUAf = fvc::interpolate(rUA);

    tmp<fvScalarMatrix> pdEqnComp;

    if (transonic)
    {
        pdEqnComp =
            (fvm::ddt(pd) + fvm::div(phi, pd) - fvm::Sp(fvc::div(phi), pd));
    }
    else
    {
        pdEqnComp =
            (fvm::ddt(pd) + fvc::div(phi, pd) - fvc::Sp(fvc::div(phi), pd));
    }


    U = rUA*UEqn.H();

    surfaceScalarField phiU
    (
        "phiU",
        (fvc::interpolate(U) & mesh.Sf()) + fvc::ddtPhiCorr(rUA, rho, U, phi)
    );

    phi = phiU +
        (
            fvc::interpolate(interface.sigmaK())*fvc::snGrad(alpha1)
          - ghf*fvc::snGrad(rho)
        )*rUAf*mesh.magSf();

    for(int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        fvScalarMatrix pdEqnIncomp
        (
            fvc::div(phi)
          - fvm::laplacian(rUAf, pd)
        );

        solve
        (
            (
                max(alpha1, scalar(0))*(psi1/rho1)
              + max(alpha2, scalar(0))*(psi2/rho2)
            )
           *pdEqnComp()
          + pdEqnIncomp
        );

        if (nonOrth == nNonOrthCorr)
        {
            dgdt =
                (pos(alpha2)*(psi2/rho2) - pos(alpha1)*(psi1/rho1))
               *(pdEqnComp & pd);
            phi += pdEqnIncomp.flux();
        }
    }

    U += rUA*fvc::reconstruct((phi - phiU)/rUAf);
    U.correctBoundaryConditions();

    p = max
        (
            (pd + gh*(alpha1*rho10 + alpha2*rho20))/(1.0 - gh*(alpha1*psi1 + alpha2*psi2)),
            pMin
        );

    rho1 = rho10 + psi1*p;
    rho2 = rho20 + psi2*p;

    Info<< "max(U) " << max(mag(U)).value() << endl;
    Info<< "min(pd) " << min(pd).value() << endl;
}
