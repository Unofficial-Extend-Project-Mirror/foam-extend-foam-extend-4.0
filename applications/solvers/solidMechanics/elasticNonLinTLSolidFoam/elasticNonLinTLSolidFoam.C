/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2007 Hrvoje Jasak
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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Application
    elasticNonLinTLSolidFoam

Description
    Finite volume structural solver employing a total strain total
    Lagrangian approach.

    Valid for finite strains, finite displacements and finite rotations.

Author
    Micheal Leonard
    Philip Cardiff

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "constitutiveModel.H"
#include "solidInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
# include "setRootCase.H"
# include "createTime.H"
# include "createMesh.H"
# include "createFields.H"
# include "createSolidInterfaceNonLin.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

  Info<< "\nStarting time loop\n" << endl;

  while(runTime.loop())
    {
      Info<< "Time: " << runTime.timeName() << nl << endl;

#     include "readSolidMechanicsControls.H"

      int iCorr = 0;
      scalar initialResidual = 0;
      lduMatrix::solverPerformance solverPerf;
      scalar relativeResidual = 1.0;
      lduMatrix::debug = 0;

      do
      {
          U.storePrevIter();

          surfaceTensorField shearGradU =
              ((I - n*n)&fvc::interpolate(gradU));

          fvVectorMatrix UEqn
              (
                  fvm::d2dt2(rho, U)
                  ==
                  fvm::laplacian(2*muf + lambdaf, U, "laplacian(DU,U)")
                  // + fvc::div(
                  //         -( (mu + lambda) * gradU )
                  //         + ( mu * gradU.T() )
                  //         + ( mu * (gradU & gradU.T()) )
                  //         + ( lambda * tr(gradU) * I )
                  //         + ( 0.5 * lambda * tr(gradU & gradU.T()) * I )
                  //         + ( sigma & gradU ),
                  //         "div(sigma)"
                  //         )
                  + fvc::div(
                      mesh.magSf()
                      *(
                          - (muf + lambdaf)*(fvc::snGrad(U)&(I - n*n))
                          + lambdaf*tr(shearGradU&(I - n*n))*n
                          + muf*(shearGradU&n)
                          + muf * (n & fvc::interpolate(gradU & gradU.T()))
                          + 0.5*lambdaf
                          *(n * tr(fvc::interpolate(gradU & gradU.T())))
                          + (n & fvc::interpolate( sigma & gradU ))
                          )
                      )
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

          //gradU = solidInterfacePtr->grad(U);
          gradU = fvc::grad(U);

#         include "calculateEpsilonSigma.H"
#         include "calculateRelativeResidual.H"

          Info<< "\tTime " << runTime.value()
              << ", Corrector " << iCorr
              << ", Solving for " << U.name()
              << " using " << solverPerf.solverName()
              << ", residual = " << solverPerf.initialResidual()
              << ", relative residual = " << relativeResidual
              << ", inner iterations " << solverPerf.nIterations() << endl;
      }
      while
          (
              solverPerf.initialResidual() > convergenceTolerance
              //relativeResidual > convergenceTolerance
              &&
              ++iCorr < nCorr
              );

      Info<< nl << "Time " << runTime.value() << ", Solving for " << U.name()
          << ", Initial residual = " << initialResidual
          << ", Final residual = " << solverPerf.initialResidual()
          << ", Relative residual = " << relativeResidual
          << ", No outer iterations " << iCorr
          << nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
          << "  ClockTime = " << runTime.elapsedClockTime() << " s"
          << endl;

#     include "writeFields.H"

      Info<< "ExecutionTime = "
          << runTime.elapsedCpuTime()
          << " s\n\n" << endl;
    }

  Info<< "End\n" << endl;

  return(0);
}


// ************************************************************************* //
