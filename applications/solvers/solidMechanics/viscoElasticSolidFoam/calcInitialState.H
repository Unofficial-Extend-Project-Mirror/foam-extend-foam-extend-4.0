
    {
        Info<< "Time: " << runTime.timeName() << nl << endl;

#       include "readSolidMechanicsControls.H"

        volScalarField mu = rheology.mu(0);
        volScalarField lambda = rheology.lambda(0);

        Info << "mu = " << average(mu.internalField()) << endl;
        Info << "lambda = " << average(lambda.internalField()) << endl;

        int iCorr = 0;
        lduSolverPerformance solverPerf;
        scalar initialResidual = 0;
        scalar err = GREAT;

        do
        {
            DU.storePrevIter();

            fvVectorMatrix DUEqn
            (
                fvm::laplacian(2*mu+lambda, DU, "laplacian(DDU,DU)")
              + fvc::div
                (
                    mu*gradDU.T()
                  + lambda*(I*tr(gradDU))
                  - (mu + lambda)*gradDU,
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
        }
        while
        (
            solverPerf.initialResidual() > convergenceTolerance
         && ++iCorr < nCorr
        );

        Info << "Solving for " << DU.name() << " using "
            << solverPerf.solverName() << " solver"
            << ", Initial residula = " << initialResidual
            << ", Final residual = " << solverPerf.initialResidual()
            << ", No outer iterations " << iCorr
            << ", Relative error: " << err << endl;

        U += DU;

#       include "calculateDSigma.H"

        sigma += DSigma;

        {
            DSigmaCorr =
                dimensionedSymmTensor
                (
                    "zero",
                    dimForce/dimArea,
                    symmTensor::zero
                );
            scalar t = runTime.value();
            scalar tNext = t + runTime.deltaT().value();
            DSigmaCorr += 2.0*rheology.mu(tNext)*Depsilon
                + rheology.lambda(tNext)*(I*tr(Depsilon));
            DSigmaCorr -= sigma;
        }

#       include "calculateStress.H"

#       include "writeHistory.H"

        Depsilon.write();

        Info<< "ExecutionTime = "
            << runTime.elapsedCpuTime()
            << " s\n\n" << endl;
    }
