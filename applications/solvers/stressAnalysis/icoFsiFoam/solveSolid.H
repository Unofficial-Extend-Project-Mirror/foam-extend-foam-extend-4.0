{
#   include "readStressedFoamControls.H"

    int iCorr = 0;
    scalar initialResidual = 0;

    do
    {
        volTensorField gradU = fvc::grad(Usolid);

        fvVectorMatrix UEqn
        (
            fvm::d2dt2(Usolid)
         ==
            fvm::laplacian(2*mu + lambda, Usolid, "laplacian(DU,U)")
          + fvc::div
            (
                mu*gradU.T() + lambda*(I*tr(gradU)) - (mu + lambda)*gradU,
                "div(sigma)"
            )
        );

        initialResidual = UEqn.solve().initialResidual();

    } while (initialResidual > convergenceTolerance && ++iCorr < nCorr);
}
