/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Application
    elasticAcpSolidFoam

Description
    Arbitrary crack propagation solver.
    Cracks may propagate along any mesh internal face.

    Please cite:
    Carolan D, Tuković Z, Murphy N, Ivankovic A, Arbitrary crack propagation
    in multi-phase materials using the finite volume method, Computational
    Materials Science, 2013, http://dx.doi.org/10.1016/j.commatsci.2012.11.049.

Author
    Zeljko Tukovic, FSB Zagreb
    Declan Carolan UCD
    Philip Cardiff UCD

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "constitutiveModel.H"
//#include "componentReferenceList.H"
#include "crackerFvMesh.H"
#include "processorPolyPatch.H"
#include "SortableList.H"
#include "solidInterface.H"
#include "solidCohesiveFvPatchVectorField.H"
#include "solidCohesiveFixedModeMixFvPatchVectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createCrackerMesh.H"
#   include "createFields.H"
#   include "createCrack.H"
//#   include "createReference.H"
#   include "createHistory.H"
#   include "readDivSigmaExpMethod.H"
#   include "createSolidInterfaceNoModify.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    lduMatrix::debug = 0;

    scalar maxEffTractionFraction = 0;

    // time rates for predictor
    volTensorField gradV = fvc::ddt(gradU);
    surfaceVectorField snGradV =
        (snGradU - snGradU.oldTime())/runTime.deltaT();

    //# include "initialiseSolution.H"

    while (runTime.run())
    {
#       include "readSolidMechanicsControls.H"
#       include "setDeltaT.H"

        runTime++;

        Info << "Time = " << runTime.timeName() << nl << endl;

        volScalarField rho = rheology.rho();
        volScalarField mu = rheology.mu();
        volScalarField lambda = rheology.lambda();
        surfaceScalarField muf = fvc::interpolate(mu);
        surfaceScalarField lambdaf = fvc::interpolate(lambda);

        if (solidInterfaceCorr)
        {
            solidInterfacePtr->modifyProperties(muf, lambdaf);
        }

        //#     include "waveCourantNo.H"

        int iCorr = 0;
        lduMatrix::solverPerformance solverPerf;
        scalar initialResidual = 0;
        scalar relativeResidual = 1;
        //scalar forceResidual = 1;
        label nFacesToBreak = 0;
        label nCoupledFacesToBreak = 0;
        bool topoChange = false;

        // Predictor step using time rates
        if (predictor)
        {
            Info<< "Predicting U, gradU and snGradU using velocity"
                << endl;
            U += V*runTime.deltaT();
            gradU += gradV*runTime.deltaT();
            snGradU += snGradV*runTime.deltaT();
        }

        do
        {
            surfaceVectorField n = mesh.Sf()/mesh.magSf();
            do
            {
                U.storePrevIter();

#               include "calculateDivSigmaExp.H"

                fvVectorMatrix UEqn
                (
                    rho*fvm::d2dt2(U)
                 ==
                    fvm::laplacian(2*muf + lambdaf, U, "laplacian(DU,U)")
                  + divSigmaExp
                );

                //#             include "setReference.H"

                if(solidInterfaceCorr)
                {
                    solidInterfacePtr->correct(UEqn);
                }

                if (relaxEqn)
                {
                    UEqn.relax();
                }

                solverPerf = UEqn.solve();

                if (aitkenRelax)
                {
#                   include "aitkenRelaxation.H"
                }
                else
                {
                    U.relax();
                }

                if (iCorr == 0)
                {
                    initialResidual = solverPerf.initialResidual();
                    aitkenInitialRes = gMax(mag(U.internalField()));
                }

                //gradU = solidInterfacePtr->grad(U);
                // use leastSquaresSolidInterface grad scheme
                gradU = fvc::grad(U);

#               include "calculateRelativeResidual.H"
                //#             include "calculateForceResidual.H"

                if (iCorr % infoFrequency == 0)
                {
                    Info<< "\tTime " << runTime.value()
                        << ", Corr " << iCorr
                        << ", Solving for " << U.name()
                        << " using " << solverPerf.solverName()
                        << ", res = " << solverPerf.initialResidual()
                        << ", rel res = " << relativeResidual;
                    if (aitkenRelax)
                    {
                        Info << ", aitken = " << aitkenTheta;
                    }
                    Info << ", inner iters " << solverPerf.nIterations() << endl;
                }
            }
            while
            (
                //iCorr++ == 0
                iCorr++ < 2
                ||
                (
                    solverPerf.initialResidual() > convergenceTolerance
                    //relativeResidual > convergenceTolerance
                 && iCorr < nCorr
                )
            );

            Info<< "Solving for " << U.name() << " using "
                << solverPerf.solverName()
                << ", Initial residual = " << initialResidual
                << ", Final residual = " << solverPerf.initialResidual()
                << ", No outer iterations " << iCorr
                << ", Relative residual " << relativeResidual << endl;

#           include "calculateTraction.H"
#           include "updateCrack.H"

            Info<< "Max effective traction fraction: "
                << maxEffTractionFraction << endl;

            // reset counter if faces want to crack
            if ((nFacesToBreak > 0)  || (nCoupledFacesToBreak > 0)) iCorr = 0;
        }
        while( (nFacesToBreak > 0)  || (nCoupledFacesToBreak > 0));

        if (cohesivePatchUPtr)
        {
            if (returnReduce(cohesivePatchUPtr->size(), sumOp<label>()))
            {
                cohesivePatchUPtr->cracking();
            }
        }
        else
        {
            if
            (
                returnReduce
                (
                    cohesivePatchUFixedModePtr->size(),
                    sumOp<label>()
                )
            )
            {
                Pout << "Number of faces in crack: "
                    << cohesivePatchUFixedModePtr->size() << endl;
                cohesivePatchUFixedModePtr->relativeSeparationDistance();
            }
        }

        // update time rates for predictor
        if (predictor)
        {
            V = fvc::ddt(U);
            gradV = fvc::ddt(gradU);
            snGradV = (snGradU - snGradU.oldTime())/runTime.deltaT();
        }

#       include "calculateEpsilonSigma.H"
#       include "writeFields.H"
#       include "writeHistory.H"

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s\n\n"
            << endl;
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
