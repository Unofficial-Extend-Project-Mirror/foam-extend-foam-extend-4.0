    {
#       include "readStressedFoamControls.H"

        int iCorr = 0;
        lduMatrix::solverPerformance solverPerf;
        scalar initialResidual = 0;

        lduMatrix::debug = 0;

#       include "backwardCoeffs.H"

        do
        {
            DU.storePrevIter();

            fvVectorMatrix DUEqn
            (
                Cn*rho*fvm::ddt(DU)
              - Co*rho*DV.oldTime()
              + Coo*rho*DV.oldTime().oldTime()
             ==
                fvm::laplacian(2*mu + lambda, DU, "laplacian(DDU,DU)")
              - fvc::laplacian(mu + lambda, DU, "laplacian(DDU,DU)")
              + fvc::div
                (
                    mu*gradDU.T()
                  + lambda*(I*tr(gradDU))
                  + mu*(gradDU&gradDU.T())
                  + 0.5*lambda*(I*tr(gradDU & gradDU.T()))
                  + (sigma & DF.T())
                  + (DSigma & DF.T()),
                    "div(sigma)"
                )
            );

            solverPerf = DUEqn.solve();

            DU.relax();

            if(iCorr == 0)
            {
                initialResidual = solverPerf.initialResidual();
            }

            gradDU = fvc::grad(DU);

            DF = gradDU.T();

#           include "calculateDSigma.H"
        } 
        while
        (
            solverPerf.initialResidual() > convergenceTolerance 
         && ++iCorr < nCorr
        );

        Info << "Solving for " << DU.name() 
            << ", Initial residual = " << initialResidual 
            << ", Final residual = " << solverPerf.initialResidual()
            << ", No outer iterations " << iCorr << endl;

        DV = fvc::ddt(DU);

        lduMatrix::debug = 1;
    }

