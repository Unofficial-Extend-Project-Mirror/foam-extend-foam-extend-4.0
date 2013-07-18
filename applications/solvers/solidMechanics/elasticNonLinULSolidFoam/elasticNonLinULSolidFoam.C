/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    elasticNonLinULSolidFoam

Description
    Finite volume structural solver employing a incremental strain updated
    Lagrangian approach.

    Valid for small strains, finite displacements and finite rotations.

Author
    Philip Cardiff


\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "rheologyModel.H"
#include "solidInterface.H"
#include "volPointInterpolation.H"
#include "pointPatchInterpolation.H"
#include "primitivePatchInterpolation.H"
#include "pointFields.H"
#include "plane.H"
#include "meshSearch.H"
#include "twoDPointCorrector.H"
#include "leastSquaresVolPointInterpolation.H"
#include "processorFvPatchFields.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"

#   include "createTime.H"

#   include "createMesh.H"

#   include "createFields.H"

#   include "readDivDSigmaExpMethod.H"

#   include "readDivDSigmaLargeStrainExpMethod.H"

#   include "readMoveMeshMethod.H"

#   include "createSolidInterface.H"

//* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info << "\nStarting time loop\n" << endl;

    for (runTime++; !runTime.end(); runTime++)
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

#       include "readStressedFoamControls.H"

        int iCorr = 0;
        lduMatrix::solverPerformance solverPerf;
        scalar initialResidual = 0;
        scalar relativeResidual = GREAT;
        lduMatrix::debug = 0;

        do
        {
            DU.storePrevIter();

            divDSigmaLargeStrainExp.storePrevIter();

#           include "calculateDivDSigmaExp.H"

#           include "calculateDivDSigmaLargeStrainExp.H"

            //----------------------------------------------------//
            //- updated lagrangian large strain momentum equation
            //----------------------------------------------------//
            fvVectorMatrix DUEqn
            (
                fvm::d2dt2(rho,DU)
              ==
                fvm::laplacian(2*muf + lambdaf, DU, "laplacian(DDU,DU)")
              + divDSigmaExp
              + divDSigmaLargeStrainExp
         );

         if(solidInterfaceCorr)
         {
             solidInterfacePtr->correct(DUEqn);
         }

         solverPerf = DUEqn.solve();

         if(iCorr == 0)
         {
             initialResidual = solverPerf.initialResidual();
         }

         DU.relax();

         if(solidInterfaceCorr)
         {
             gradDU = solidInterfacePtr->grad(DU);
         }
         else
         {
             gradDU = fvc::grad(DU);
         }

#        include "calculateDEpsilonDSigma.H"

#        include "calculateRelativeResidual.H"

         Info << "\tTime " << runTime.value()
             << ", Corrector " << iCorr
             << ", Solving for " << DU.name()
             << " using " << solverPerf.solverName()
             << ", residual = " << solverPerf.initialResidual()
             << ", residualDU = " << relativeResidual
             << ", inner iterations " << solverPerf.nIterations() << endl;
      }
      while
      (
          //solverPerf.initialResidual() > convergenceTolerance
          relativeResidual > convergenceTolerance
          && ++iCorr < nCorr
      );

     lduMatrix::debug = 1;

     Info << nl << "Time " << runTime.value() << ", Solving for " << DU.name()
         << ", Initial residual = " << initialResidual
         << ", Final residual = " << solverPerf.initialResidual()
         << ", No outer iterations " << iCorr << endl;

#    include "rotateFields.H"

#    include "moveMesh.H"

#    include "writeFields.H"

     //- total force
     forAll(mesh.boundary(), patchi)
     {
          vector force = sum(mesh.Sf().boundaryField()[patchi] & sigma.boundaryField()[patchi]);
          Info << "force on " << mesh.boundary()[patchi].name()
              << " is " << force << endl;
     }

     Info << nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
         << "  ClockTime = " << runTime.elapsedClockTime() << " s"
         << endl;
    }

    Info<< "End\n" << endl;

    return(0);
}

// ************************************************************************* //
