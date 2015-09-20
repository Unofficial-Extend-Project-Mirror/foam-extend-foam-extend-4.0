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
    elasticPlasticSolidFoam

Description
    Transient/steady-state segregated finite-volume solver for small strain
    elastic plastic solid bodies.

    Displacement increment field DU is solved for using a total Lagrangian
    approach, also generating the strain tensor field epsilon, the plastic
    strain field epsilonP and stress tensor field sigma.

Author
    A. Karac UCD/Zenica
    P. Cardiff UCD
    Aitken relaxation by T. Tang DTU

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "constitutiveModel.H"
#include "solidContactFvPatchVectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
#   include "createFields.H"
#   include "createHistory.H"
#   include "readDivDSigmaExpMethod.H"

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
        scalar plasticResidual = 1.0;
        lduMatrix::debug = 0;

        do
        {
            DU.storePrevIter();

#           include "calculateDivDSigmaExp.H"

            fvVectorMatrix DUEqn
            (
                rho*fvm::d2dt2(DU)
             ==
                fvm::laplacian(2*muf + lambdaf, DU, "laplacian(DDU,DU)")
              + divDSigmaExp
              - fvc::div(2*muf*(mesh.Sf() & fvc::interpolate(DEpsilonP)))
            );

            solverPerf = DUEqn.solve();

            if (iCorr == 0)
            {
                initialResidual = solverPerf.initialResidual();
            }

            if (aitkenRelax)
            {
#               include "aitkenRelaxation.H"
            }
            else
            {
                DU.relax();
            }

            gradDU = fvc::grad(DU);

#           include "calculateRelativeResidual.H"
#           include "calculateDEpsilonDSigma.H"

            // correct plastic strain increment
            rheology.correct();

#           include "calculatePlasticResidual.H"

            if (iCorr % infoFrequency == 0)
            {
                Info<< "\tTime " << runTime.value()
                    << ", Corr " << iCorr
                    //<< ", Solving for " << DU.name()
                    // << " using " << solverPerf.solverName()
                    << ", res = " << solverPerf.initialResidual()
                    << ", rel res = " << relativeResidual
                    << ", plastic res = " << plasticResidual;
                if (aitkenRelax)
                {
                    Info<< ", aitken = " << aitkenTheta;
                }
                Info<< ", inner iters = " << solverPerf.nIterations() << endl;
            }
        }
        while
        (
            iCorr++ < 2
            ||
            (
                //solverPerf.initialResidual() > convergenceTolerance
                relativeResidual > convergenceTolerance
             && iCorr < nCorr
            )
        );

        Info<< nl << "Time " << runTime.value() << ", Solving for " << DU.name()
            << ", Initial residual = " << initialResidual
            << ", Final residual = " << solverPerf.initialResidual()
            << ", Final rel residual = " << relativeResidual
            << ", No outer iterations " << iCorr << endl;

        // Update total quantities
        U += DU;
        epsilon += DEpsilon;
        epsilonP += rheology.DEpsilonP();
        sigma += DSigma;

        // Update yields stresses
        rheology.updateYieldStress();

#       include "writeFields.H"
#       include "writeHistory.H"

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

  Info<< "End\n" << endl;

  return(0);
}


// ************************************************************************* //
