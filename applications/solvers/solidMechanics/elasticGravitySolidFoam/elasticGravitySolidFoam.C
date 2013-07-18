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
    elasticGravitySolidFoam

Description
    Transient/steady-state segregated finite-volume solver for small strain
    elastic solid bodies.

    Displacement field U is solved for using a total Lagrangian approach,
    also generating the strain tensor field epsilon and stress tensor
    field sigma.

    With optional multi-material solid interface correction ensuring
    correct tractions on multi-material interfaces.

    Graivty added as a body force.

Author
    Philip Cardiff
    multi-material by Tukovic et al. 2012

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "rheologyModel.H"
#include "solidInterface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
# include "setRootCase.H"

# include "createTime.H"

# include "createMesh.H"

# include "createFields.H"

# include "readDivSigmaExpMethod.H"

# include "createSolidInterface.H"

# include "readGravity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

  Info<< "\nCalculating displacement field\n" << endl;

  while(runTime.loop())
    {
      Info<< "Time: " << runTime.timeName() << nl << endl;

#     include "readStressedFoamControls.H"

      int iCorr = 0;
      scalar initialResidual = 0;
      lduMatrix::solverPerformance solverPerf;
      scalar relativeResidual = GREAT;

      lduMatrix::debug=0;

      do
        {
	  U.storePrevIter();

#         include "calculateDivSigmaExp.H"

	  //- linear momentum equation
	  fvVectorMatrix UEqn
            (
	     fvm::d2dt2(rho, U)
	     ==
	     fvm::laplacian(2*muf + lambdaf, U, "laplacian(DU,U)")
	     + divSigmaExp
	     + rho*gravity
	     );

	  if(solidInterfaceCorr)
	    {
	      solidInterfacePtr->correct(UEqn);
	    }

	  solverPerf = UEqn.solve();

	  if(iCorr == 0)
	    {
	      initialResidual = solverPerf.initialResidual();
	    }

	  U.relax();

	  if(solidInterfaceCorr)
	    {
	      gradU = solidInterfacePtr->grad(U);
	    }
	  else
	    {
	      gradU = fvc::grad(U);
	    }

#         include "calculateRelativeResidual.H"

	  Info << "\tTime " << runTime.value()
	       << ", Corrector " << iCorr
	       << ", Solving for " << U.name()
	       << " using " << solverPerf.solverName()
	       << ", residual = " << solverPerf.initialResidual()
	       << ", relative residual = " << relativeResidual << endl;
        }
	while
	  (
	   //solverPerf.initialResidual() > convergenceTolerance
	   relativeResidual > convergenceTolerance
	   &&
	   ++iCorr < nCorr
	   );

      Info << nl << "Time " << runTime.value() << ", Solving for " << U.name()
	   << ", Initial residual = " << initialResidual
	   << ", Final residual = " << solverPerf.initialResidual()
	   << ", Relative residual = " << relativeResidual
	   << ", No outer iterations " << iCorr
	   << nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
	   << "  ClockTime = " << runTime.elapsedClockTime() << " s"
	   << endl;

      lduMatrix::debug=0;

#     include "calculateEpsilonSigma.H"

#     include "writeFields.H"

      Info<< "ExecutionTime = "
	  << runTime.elapsedCpuTime()
	  << " s\n\n" << endl;
    }

  Info<< "End\n" << endl;

  return(0);
}


// ************************************************************************* //
