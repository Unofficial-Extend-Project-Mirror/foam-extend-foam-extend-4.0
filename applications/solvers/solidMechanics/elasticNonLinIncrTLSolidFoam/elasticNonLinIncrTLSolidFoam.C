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
    elasticNonLinIncrTLSolidFoam

Description
    Finite volume structural solver employing an incremental strain total
    Lagrangian approach.
    For elastic solids.
    Valid for finite strains, finite displacements and finite rotations.

Author
    Philip Cardiff UCD
    Micheal Leonard UCD

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "constitutiveModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
#   include "createFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while(runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

#       include "readSolidMechanicsControls.H"

        int iCorr = 0;
        scalar initialResidual = 0;
        lduMatrix::solverPerformance solverPerf;
        scalar relativeResidual = 1.0;
        lduMatrix::debug=0;

        do
        {
            DU.storePrevIter();

            fvVectorMatrix DUEqn
            (
                fvm::d2dt2(rho, DU)
             ==
                fvm::laplacian(2*mu + lambda, DU, "laplacian(DDU,DU)")
              + fvc::div
                (
                  - ( (mu + lambda) * gradDU )
                  + (
                        mu *
                        (
                            gradDU.T()
                          + (gradDU & gradU.T())
                          + (gradU & gradDU.T())
                          + (gradDU & gradDU.T())
                        )
                    )
                  + ( lambda * tr(DEpsilon) * I )
                  + ( DSigma & gradU )
                  + ( (sigma + DSigma) & gradDU ),
                    "div(sigma)"
                )
            );

            solverPerf = DUEqn.solve();

            if (iCorr == 0)
            {
                initialResidual = solverPerf.initialResidual();
            }

            DU.relax();

            gradDU = fvc::grad(DU);

#           include "calculateDEpsilonDSigma.H"
#           include "calculateRelativeResidual.H"

            Info<< "\tTime " << runTime.value()
                << ", Corrector " << iCorr
                << ", Solving for " << DU.name()
                << " using " << solverPerf.solverName()
                << ", residual = " << solverPerf.initialResidual()
                << ", relative residual = " << relativeResidual
                << ", inner iterations = " << solverPerf.nIterations() << endl;
        }
        while
        (
            solverPerf.initialResidual() > convergenceTolerance
            //relativeResidual > convergenceTolerance
         && ++iCorr < nCorr
        );

        Info<< nl << "Time " << runTime.value() << ", Solving for " << DU.name()
            << ", Initial residual = " << initialResidual
            << ", Final residual = " << solverPerf.initialResidual()
            << ", Relative residual = " << relativeResidual
            << ", No outer iterations " << iCorr
            << nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << endl;

        U += DU;
        gradU += gradDU;
        epsilon += DEpsilon;
        sigma += DSigma;

#       include "writeFields.H"

        Info<< "ExecutionTime = "
            << runTime.elapsedCpuTime()
            << " s\n\n" << endl;
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
