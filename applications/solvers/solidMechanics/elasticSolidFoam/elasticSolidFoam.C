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
    elasticSolidFoam

Description
    Transient/steady-state segregated finite-volume solver for small strain
    elastic solid bodies.

    Displacement field U is solved for using a total Lagrangian approach,
    also generating the strain tensor field epsilon and stress tensor
    field sigma.

    With optional multi-material solid interface correction ensuring
    correct tractions on multi-material interfaces

Author
    Philip Cardiff
    multi-material by Tukovic et al. 2012

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "constitutiveModel.H"
#include "solidInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
#   include "createFields.H"
#   include "createHistory.H"
#   include "readDivSigmaExpMethod.H"
#   include "createSolidInterfaceNoModify.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

  while(runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

#       include "readSolidMechanicsControls.H"

        int iCorr = 0;
        lduMatrix::solverPerformance solverPerf;
        scalar initialResidual = 1.0;
        scalar relativeResidual = 1.0;
        lduMatrix::debug = 0;

        if (predictor)
        {
            Info<< "\nPredicting U, gradU and snGradU based on V,"
                << "gradV and snGradV\n" << endl;
            U += V*runTime.deltaT();
            gradU += gradV*runTime.deltaT();
            snGradU += snGradV*runTime.deltaT();
        }

        do
        {
            U.storePrevIter();

#           include "calculateDivSigmaExp.H"

            // linear momentum equation
            fvVectorMatrix UEqn
            (
                rho*fvm::d2dt2(U)
             ==
                fvm::laplacian(2*muf + lambdaf, U, "laplacian(DU,U)")
              + divSigmaExp
            );

            if (solidInterfaceCorr)
            {
                solidInterfacePtr->correct(UEqn);
            }

            solverPerf = UEqn.solve();

            if (iCorr == 0)
            {
                initialResidual = solverPerf.initialResidual();
                aitkenInitialRes = gMax(mag(U.internalField()));
            }

            if (aitkenRelax)
            {
#               include "aitkenRelaxation.H"
            }
            else
            {
                U.relax();
            }

            gradU = fvc::grad(U);

#           include "calculateRelativeResidual.H"

            if (iCorr % infoFrequency == 0)
            {
                Info<< "\tTime " << runTime.value()
                    << ", Corrector " << iCorr
                    << ", Solving for " << U.name()
                    << " using " << solverPerf.solverName()
                    << ", res = " << solverPerf.initialResidual()
                    << ", rel res = " << relativeResidual;

                if (aitkenRelax)
                {
                    Info<< ", aitken = " << aitkenTheta;
                }
                Info<< ", inner iters = " << solverPerf.nIterations() << endl;
            }
        }
        while
        (
            iCorr++ == 0
            ||
            (
                solverPerf.initialResidual() > convergenceTolerance
                //relativeResidual > convergenceTolerance
             && iCorr < nCorr
            )
        );

        Info<< nl << "Time " << runTime.value() << ", Solving for " << U.name()
            << ", Initial residual = " << initialResidual
            << ", Final residual = " << solverPerf.initialResidual()
            << ", Relative residual = " << relativeResidual
            << ", No outer iterations " << iCorr
            << nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << endl;

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
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << endl << endl;
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
