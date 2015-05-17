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
    elasticPlasticNonLinTLSolidFoam

Description
    Finite volume structural solver employing an incremental strain total
    Lagrangian approach, with Mises plasticity.

    Valid for finite strains, finite displacements and finite rotations.

Author
    Philip Cardiff UCD
    Aitken relaxation by T. Tang DTU

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "constitutiveModel.H"
#include "solidContactFvPatchVectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
# include "setRootCase.H"
# include "createTime.H"
# include "createMesh.H"
# include "createFields.H"
# include "createHistory.H"
# include "readDivDSigmaExpMethod.H"
# include "readDivDSigmaNonLinExpMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while(runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

#       include "readSolidMechanicsControls.H"

        int iCorr = 0;
        scalar initialResidual = 0;
        lduMatrix::solverPerformance solverPerf;
        scalar relativeResidual = GREAT;
        lduMatrix::debug = 0;

        do
        {
            DU.storePrevIter();

#           include "calculateDivDSigmaExp.H"
#           include "calculateDivDSigmaNonLinExp.H"

            // Incremental form of the
            // linear momentum conservation
            // ensuring conservation of total momentum
            fvVectorMatrix DUEqn
            (
                fvm::d2dt2(rho, DU)
             ==
                fvm::laplacian(2*muf + lambdaf, DU, "laplacian(DDU,DU)")
              + divDSigmaExp
              + divDSigmaNonLinExp
            //- fvc::div(2*mu*DEpsilonP, "div(sigma)")
              - fvc::div
                (
                    2*muf*(mesh.Sf() & fvc::interpolate(DEpsilonP))
                )
            );

            if (largeStrainOverRelax)
            {
                // the terms (gradDU & gradU.T()) and (gradU & gradDU.T())
                // are linearly dependent of DU and represent initial
                // displacement effect
                // which can cause convergence difficulties when treated
                // explicitly
                // so we implicitly over-relax with gradU & gradDU here
                // which tends to help convergence
                // this should improve convergence when gradU is large
                // but maybe not execution time
                DUEqn -=
                    fvm::laplacian
                    (
                        (2*mu + lambda)*gradU, DU, "laplacian(DDU,DU)"
                    )
                  - fvc::div((2*mu + lambda)*(gradU & gradDU), "div(sigma)");
            }

            if (nonLinearSemiImplicit)
            {
                // experimental
                // we can treat the nonlinear term (gradDU & gradDU.T()) in a
                // semi-implicit over-relaxed manner
                // this should improve convergence when gradDU is large
                // but maybe not execution time
                DUEqn -=
                    fvm::laplacian
                    (
                        (2*mu + lambda)*gradDU, DU, "laplacian(DDU,DU)"
                    )
                  - fvc::div((2*mu + lambda)*(gradDU & gradDU), "div(sigma)");
            }

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

            // correct plasticty term
            rheology.correct();

#           include "calculateDEpsilonDSigma.H"
#           include "calculateRelativeResidual.H"

            if(iCorr % infoFrequency == 0)
            {
                Info<< "\tTime " << runTime.value()
                    << ", Corrector " << iCorr
                    << ", Solving for " << DU.name()
                    << " using " << solverPerf.solverName()
                    << ", res = " << solverPerf.initialResidual()
                    << ", rel res = " << relativeResidual;

                if(aitkenRelax)
                {
                    Info<< ", aitken = " << aitkenTheta;
                }
                Info<< ", iters = " << solverPerf.nIterations() << endl;
            }
        }
        while
        (
            iCorr++ == 0
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
            << ", Relative residual = " << relativeResidual
            << ", No outer iterations " << iCorr
            << nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << endl;

        // Update total quantities
        U += DU;
        gradU += gradDU;
        epsilon += DEpsilon;
        epsilonP += rheology.DEpsilonP();
        sigma += DSigma;
        rheology.updateYieldStress();
        rho = rho/det(I+gradU);

#       include "writeFields.H"
#       include "writeHistory.H"

        Info<< nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << endl;
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
