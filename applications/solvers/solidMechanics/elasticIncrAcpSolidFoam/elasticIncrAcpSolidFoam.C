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
    elasticIncrAcpSolidFoam

Description
    Incremental form of elasticAcpSolidFoam
    arbitrary crack propagation solver

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
#   include "readDivDSigmaExpMethod.H"
#   include "createSolidInterfaceIncrNoModify.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    lduMatrix::debug = 0;

    scalar maxEffTractionFraction = 0;

    while (runTime.run())
    {
#       include "readSolidMechanicsControls.H"
#       include "setDeltaT.H"

        runTime++;

        Info<< "\nTime = " << runTime.timeName() << " s\n" << endl;

        volScalarField rho = rheology.rho();
        volScalarField mu = rheology.mu();
        volScalarField lambda = rheology.lambda();
        surfaceScalarField muf = fvc::interpolate(mu);
        surfaceScalarField lambdaf = fvc::interpolate(lambda);

        solidInterfacePtr->modifyProperties(muf, lambdaf);
        //#     include "waveCourantNo.H"

        int iCorr = 0;
        lduMatrix::solverPerformance solverPerf;
        scalar initialResidual = 0;
        scalar relativeResidual = 1;
        //scalar forceResidual = 1;
        label nFacesToBreak = 0;
        label nCoupledFacesToBreak = 0;
        bool topoChange = false;

        // DU from the previous timestep is usually a good guess
        // for the next timestep, but it can cause faces to prematurely
        // crack.
        // so I will reduce DU here to stop this happening
        if (!predictor)
        {
            DU *= 0.0;
        }

        do
        {
            surfaceVectorField n = mesh.Sf()/mesh.magSf();
            do
            {
                DU.storePrevIter();

#               include "calculateDivDSigmaExp.H"

                fvVectorMatrix DUEqn
                (
                    rho*fvm::d2dt2(DU)
                    ==
                    fvm::laplacian(2*muf + lambdaf, DU, "laplacian(DDU,DU)")
                  + divDSigmaExp
                );

                //#             include "setReference.H"

                if(solidInterfacePtr)
                {
                    solidInterfacePtr->correct(DUEqn);
                }

                //DUEqn.relax();

                solverPerf = DUEqn.solve();

                if (aitkenRelax)
                {
#                   include "aitkenRelaxation.H"
                }
                else
                {
                    DU.relax();
                }

                if (iCorr == 0)
                {
                    initialResidual = solverPerf.initialResidual();
                    aitkenInitialRes = gMax(mag(DU.internalField()));
                }

                //gradDU = solidInterfacePtr->grad(DU);
                // use leastSquaresSolidInterface grad scheme
                gradDU = fvc::grad(DU);

#               include "calculateRelativeResidual.H"

                if (iCorr % infoFrequency == 0)
                {
                    Info<< "\tTime " << runTime.value()
                        << ", Corr " << iCorr
                        << ", Solving for " << DU.name()
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
                iCorr++ < 10
                ||
                (
                    //solverPerf.initialResidual() > convergenceTolerance
                    relativeResidual > convergenceTolerance
                    &&
                    iCorr < nCorr
                )
            );

            Info<< "Solving for " << DU.name() << " using "
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

        if (cohesivePatchDUPtr)
        {
            if (returnReduce(cohesivePatchDUPtr->size(), sumOp<label>()))
            {
                cohesivePatchDUPtr->cracking();
            }
        }
        else
        {
            if
            (
                returnReduce
                (
                    cohesivePatchDUFixedModePtr->size(),
                    sumOp<label>()
                )
            )
            {
                Pout << "Number of faces in crack: "
                    << cohesivePatchDUFixedModePtr->size() << endl;
                cohesivePatchDUFixedModePtr->relativeSeparationDistance();
            }
        }

#       include "calculateDEpsilonDSigma.H"

        // update total quantities
        U += DU;
        epsilon += DEpsilon;
        sigma += DSigma;

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
