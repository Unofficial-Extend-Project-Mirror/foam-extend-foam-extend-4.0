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
    viscoElasticSolidFoam

Description
    visco-elastic small strain solver using finite volume method,
    using an incremental approach

Author
    Zeljko Tukovic FSB Zagreb

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "constitutiveModel.H"
//#include "componentReferenceList.H"
//#include "patchToPatchInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
#   include "createFields.H"
#   include "createHistory.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    Info<< "Note: the results must be written for every time-step"
        << " as they are used to calculate the current stress" << endl;

    lduMatrix::debug = 0;
    scalar m = 0.5;
    surfaceVectorField n = mesh.Sf()/mesh.magSf();

    while(runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

#       include "readSolidMechanicsControls.H"

        volScalarField mu = rheology.mu(m*runTime.deltaT().value());
        volScalarField lambda = rheology.lambda(m*runTime.deltaT().value());
        surfaceScalarField muf = fvc::interpolate(mu);
        surfaceScalarField lambdaf = fvc::interpolate(lambda);
        Info << "average mu = " << average(muf.internalField()) << endl;
        Info << "average lambda = " << average(lambdaf.internalField()) << endl;

        int iCorr = 0;
        lduMatrix::solverPerformance solverPerf;
        scalar initialResidual = 1.0;
        scalar residual = 1.0;
        surfaceSymmTensorField DSigmaCorrf = fvc::interpolate(DSigmaCorr);

        //label nCrackedFaces = 0;
        // cracking loop if you use cohesive boundaries
        //do
        //{
        do
        {
            surfaceTensorField sGradDU =
                (I - n*n)&fvc::interpolate(gradDU);

            DU.storePrevIter();

            fvVectorMatrix DUEqn
            (
                rho*fvm::d2dt2(DU)
                ==
                fvm::laplacian(2*muf+lambdaf, DU, "laplacian(DDU,DU)")
                + fvc::div
                (
                    mesh.magSf()*
                    (
                       - (muf + lambdaf)*(fvc::snGrad(DU)&(I - n*n))
                       + lambdaf*tr(sGradDU&(I - n*n))*n
                       + muf*(sGradDU&n)
                       + (n&DSigmaCorrf)
                    )
                )
            );

//            // add an increment of gravity on the first time-step
//            if (runTime.timeIndex() == 1)
//            {
//                DUEqn -= (rho*g);
//            }

            solverPerf = DUEqn.solve();

            DU.relax();

            if (iCorr == 0)
            {
                initialResidual = solverPerf.initialResidual();
            }

            gradDU = fvc::grad(DU);

#           include "calculateDSigma.H"
#           include "calcResidual.H"

            if (iCorr % infoFrequency == 0)
            {
                Info<< "\tTime " << runTime.value()
                    << ", Corrector " << iCorr
                    << ", Solving for " << U.name()
                    << " using " << solverPerf.solverName()
                    << ", res = " << solverPerf.initialResidual()
                    << ", rel res = " << residual
                    << ", inner iters = " << solverPerf.nIterations() << endl;
            }
        }
        while
        (
            // solverPerf.initialResidual() > convergenceTolerance
            residual > convergenceTolerance
         && ++iCorr < nCorr
        );

        Info<< "Solving for " << DU.name() << " using "
            << solverPerf.solverName() << " solver"
            << ", Initial residula = " << initialResidual
            << ", Final residual = " << solverPerf.initialResidual()
            << ", No outer iterations " << iCorr
            << ", Relative error: " << residual << endl;

        //#           include "updateCrack.H"
        //}
        //while(nCrackedFaces > 0);

        U += DU;

#       include "calculateSigma.H"
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
