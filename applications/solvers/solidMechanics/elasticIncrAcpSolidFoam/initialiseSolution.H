if (initialiseSolution)
{
    Info << "Initialise solution" << endl;

    const dictionary& stressControl =
        mesh.solutionDict().subDict("stressedFoam");

    int nCorr(readInt(stressControl.lookup("nCorrectors")));
    scalar convergenceTolerance
    (
        readScalar(stressControl.lookup("U"))
    );
    
    volScalarField rho = rheology.rho();
    surfaceScalarField mu = fvc::interpolate(rheology.mu());
    surfaceScalarField lambda = fvc::interpolate(rheology.lambda());
    //interface.modifyProperties(mu, lambda);
    interfacePtr->modifyProperties(mu, lambda);

    int iCorr = 0;
    lduMatrix::solverPerformance solverPerf;
    scalar initialResidual = 0;
    scalar relRes = GREAT;

    surfaceVectorField n = mesh.Sf()/mesh.magSf();

    do
    {
        U.storePrevIter();

        surfaceTensorField sGradU =
            ((I - n*n)&fvc::interpolate(gradU));

        fvVectorMatrix UEqn
        (
            fvm::laplacian(2*mu + lambda, U, "laplacian(DU,U)")
          + fvc::div
            (
                mesh.magSf()
               *(
                  - (mu + lambda)*(snGradU & (I - n*n))
                  + mu*(sGradU & n)
                  + lambda*tr(sGradU&(I - n*n))*n
                )
            )
        );

#       include "setReference.H"

	interfacePtr->correct(UEqn);

        solverPerf = UEqn.solve();

        U.relax();

        if(iCorr == 0)
        {
            initialResidual = solverPerf.initialResidual();
        }

        //gradU = fvc::grad(U);
        //gradU = interface.grad(U);
        gradU = interfacePtr->grad(U);
        snGradU = fvc::snGrad(U);

        relRes =
            gMax
            (
                mag(U.internalField() - U.prevIter().internalField())
               /(
                    mag
                    (
                        U.internalField()
                      - U.oldTime().internalField()
                    )
                  + SMALL
                )
            );
    }
    while
    (
        relRes > convergenceTolerance
     && ++iCorr < nCorr
    );

    Info << "Solving for " << U.name() << " using "
        << solverPerf.solverName()
        << ", Initial residual = " << initialResidual
        << ", Final residual = " << solverPerf.initialResidual()
        << ", No outer iterations " << iCorr
        << ", Relative residual " << relRes << endl;
    
#   include "calculateTraction.H"

    scalarField effTraction = 
      mag(n.internalField() & traction.internalField());
        // cohesivePatchU.effectiveTraction
        // (
        //     traction.internalField(),
        //     n.internalField()
        // );

    scalar maxEffTraction = 0.0;

    if
    (
        crackPropagationFromSpecifiedPatches
     && crackPropagationPatchesInternalFaces.size()
    )
    {
        maxEffTraction =
            gMax
            (
                scalarField
                (
                    effTraction, 
                    crackPropagationPatchesInternalFaces
                )
            );
    }
    else
    {
        maxEffTraction = gMax(cohesiveZone.internalField()*effTraction);
    }

    Info << "Max effTraction: " << maxEffTraction << endl;

//     if (maxEffTraction >= CzLaw.sigmaMax().value())
//     {
//         FatalErrorIn(args.executable())
//             << "Maximal effective traction is greater then critical traction" 
//                 << abort(FatalError);
//     }

    U.oldTime() = U;
}
