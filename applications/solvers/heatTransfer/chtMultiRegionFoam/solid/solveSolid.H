{
    for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
    {
        solve
        (
            fvm::ddt(rhosCps[i], Ts[i]) - fvm::laplacian(Ks[i], Ts[i])
        );
    }

    Info<< "Min/max T:" << min(Ts[i]) << ' ' << max(Ts[i])
        << endl;
}
