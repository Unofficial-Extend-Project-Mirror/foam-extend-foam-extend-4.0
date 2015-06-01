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
    elasticOrthoSolidFoam

Description
    Transient/steady-state segregated finite-volume solver for small strain
    elastic orthotropic solid bodies allowing for general principal material
    directions.

    Please cite:
    Cardiff P, Karac A & Ivankovic A, A Large Strain Finite Volume Method for
    Orthotropic Bodies with General Material Orientations, Computer Methods
    in Applied Mechanics & Engineering, 2013,
    http://dx.doi.org/10.1016/j.cma.2013.09.008

Author
    Philip Cardiff UCD

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
#   include "createSolidInterfaceOrthotropic.H"

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

        do
        {
            U.storePrevIter();

#           include "calculateDivSigmaExp.H"

            //- Linear momentum equation
            fvVectorMatrix UEqn
            (
                rho*fvm::d2dt2(U)
             ==
                fvm::laplacian(Kf, U, "laplacian(K,U)")
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
            }

            U.relax();

            gradU = fvc::grad(U); // use leastSquaresSolidInterface

            //#         include "setPlaneStressGradU.H"

#           include "calculateRelativeResidual.H"

            if (iCorr % infoFrequency == 0)
            {
                Info<< "\tTime " << runTime.value()
                    << ", Corr " << iCorr
                    << ", Solving for " << U.name()
                    << " using " << solverPerf.solverName()
                    << ", res = " << solverPerf.initialResidual()
                    << ", rel res = " << relativeResidual
                    << ", inner iters " << solverPerf.nIterations() << endl;
            }
        }
        while
        (
            solverPerf.initialResidual() > convergenceTolerance
         && ++iCorr < nCorr
        );

        Info<< nl << "Time " << runTime.value() << ", Solving for " << U.name()
            << ", Initial residual = " << initialResidual
            << ", Final residual = " << solverPerf.initialResidual()
            << ", No outer iterations " << iCorr
            << nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << endl;

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
