psim1 = psi;
psi = psi_tmp;

// Resolution of the linear system.
fvScalarMatrix psiEqn
(
    Crel*fvm::ddt(psi)
 ==
    fvm::laplacian(Krel, psi, "laplacian(Krel,psi)")
  + gradkz
);
psiEqn.relax();
psiEqn.solve();

// Update of the varying transport properties.
thtil = 0.5*
(
    (1 + sign(psi)) + (1 - sign(psi))*
    pow((1+pow(mag(alpha*psi),n)),-(1-(1/n)))
);
Krel = 0.5*
(
    (1 + sign(psi))*K + (1 - sign(psi))*K*pow(thtil, 0.5)*
    pow((1 - pow((1 - pow(thtil, (n/(n - 1)))),(1 - (1/n)))), 2)
);
Crel= S + 0.5*
(
    (1 - sign(psi))*
    (
        (thetas - thetar)*(thtil - thtil_tmp)*
        (
            1./((usf*pos(psi - psi_tmp)*pos(psi_tmp - psi))
            + psi - psi_tmp)
        )
    )
);

// Update of the gravity term.
gradk = fvc::grad(Krel);
gradkz = gradk.component(2);

// Computation of the field of residuals for the exit test of
// the Picard loop.
err = psi - psim1;

// Computation of the residual for the exit test of the Picard
// loop. Note the use of gSum instead of sum.
crit = gSumMag(err)/nbMesh;

// exit test for the Picard loop.
if (crit < precPicard)
{
    Info << " Erreur = " << crit
        << " Picard = " << pic
        << nl << endl;

    currentPicard=pic;

    break;
}
