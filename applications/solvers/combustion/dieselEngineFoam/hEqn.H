{
    solve
    (
        fvm::ddt(rho, h)
      + mvConvection->fvmDiv(phi, h)
      - fvm::laplacian(turbulence->alphaEff(), h)
     ==
       DpDt
     + dieselSpray.heatTransferSource()
    );

    thermo->correct();

    forAll(dQ, i)
    {
        dQ[i] = 0.0;
    }

    scalarField cp(dQ.size(), 0.0);

    forAll(Y, i)
    {
        volScalarField RRi = chemistry.RR(i);

        forAll(h, celli)
        {
            scalar Ti = T[celli];
            cp[celli] += Y[i][celli]*chemistry.specieThermo()[i].Cp(Ti);
            scalar hi = chemistry.specieThermo()[i].h(Ti);
            scalar RR = RRi[celli];
            dQ[celli] -= hi*RR;
        }

    }

    forAll(dQ, celli)
    {
        dQ[celli] /= cp[celli];
    }
}
