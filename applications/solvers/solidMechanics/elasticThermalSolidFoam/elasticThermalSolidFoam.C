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
    elasticThermalSolidFoam

Description
    Transient/steady-state segregated finite-volume solver for small strain
    elastic thermal solid bodies. Temperature is solved and then coupled
    displacement is solved.

    Displacement field U is solved for using a total Lagrangian approach,
    also generating the strain tensor field epsilon and stress tensor
    field sigma and temperature field T.

Author
   Philip Cardiff UCD

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "constitutiveModel.H"
#include "thermalModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
# include "setRootCase.H"
# include "createTime.H"
# include "createMesh.H"
# include "createFields.H"
# include "readDivSigmaExpMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

  Info<< "\nStarting time loop\n" << endl;

  while(runTime.loop())
    {
      Info<< "Time: " << runTime.timeName() << nl << endl;
      
#     include "readSolidMechanicsControls.H"
      
      int iCorr = 0;
      scalar initialResidual = 1.0;
      scalar relResT = 1.0;
      scalar relResU = 1.0;
      lduMatrix::solverPerformance solverPerfU;
      lduMatrix::solverPerformance solverPerfT;
      lduMatrix::debug = 0;

      // solve energy equation for temperature
      // the loop is for non-orthogonal corrections
      Info << "Solving for " << T.name() << nl;
      do
	{
	  T.storePrevIter();

	  fvScalarMatrix TEqn
	    (
	     rhoC*fvm::ddt(T) == fvm::laplacian(k, T, "laplacian(k,T)")
	     );
	  
	  solverPerfT = TEqn.solve();

	  T.relax();

#         include "calculateRelResT.H"

          if(iCorr % infoFrequency == 0)
            {
	      Info << "\tCorrector " << iCorr
		   << ", residual = " << solverPerfT.initialResidual()
		   << ", relative res = " << relResT
		   << ", inner iters = " << solverPerfT.nIterations() << endl;
	    }
	}
      while
	(
	 relResT > convergenceToleranceT
	 &&
	 ++iCorr < nCorr
	 );

      Info << "Solved for " << T.name()
	   << " using " << solverPerfT.solverName()
	   << " in " << iCorr << " iterations"
	   << ", residual = " << solverPerfT.initialResidual()
	   << ", relative res = " << relResT << nl
	   << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
	   << ", ClockTime = " << runTime.elapsedClockTime() << " s" 
	   << endl;
	
      // solve momentum equation for displacement
      iCorr = 0;
      volVectorField gradThreeKalphaDeltaT = fvc::grad(threeKalpha*(T-T0), "grad(threeKalphaDeltaT)");
      surfaceVectorField threeKalphaDeltaTf = mesh.Sf()*threeKalphaf*fvc::interpolate(T-T0, "deltaT");

      Info << "Solving for " << U.name() << nl;
      do
        {
	  U.storePrevIter();

#         include "calculateDivSigmaExp.H"

	  // linear momentum equaiton
	  fvVectorMatrix UEqn
            (
	     rho*fvm::d2dt2(U)
	     ==
	     fvm::laplacian(2*muf + lambdaf, U, "laplacian(DU,U)")
	     + divSigmaExp
	     );

	  solverPerfU = UEqn.solve();

	  if(aitkenRelax)
            {
#             include "aitkenRelaxation.H"
            }
          else
            {
              U.relax();
            }

	  gradU = fvc::grad(U);

#         include "calculateRelResU.H"

	  if(iCorr == 0)
	    {
	      initialResidual = solverPerfU.initialResidual();
	    }

          if(iCorr % infoFrequency == 0)
            {
	      Info << "\tCorrector " << iCorr
		   << ", residual = " << solverPerfU.initialResidual()
		   << ", relative res = " << relResU;
              if(aitkenRelax) Info << ", aitken = " << aitkenTheta;
              Info << ", inner iters = " << solverPerfU.nIterations() << endl;
	    }
	}
      while
	(
         iCorr++ == 0
         ||
         (//solverPerfU.initialResidual() > convergenceTolerance
	  relResU > convergenceToleranceU
	  &&
          iCorr < nCorr)
	 );
	
      Info << "Solved for " << U.name()
	   << " using " << solverPerfU.solverName()
	   << " in " << iCorr << " iterations"
	   << ", initial res = " << initialResidual
	   << ", final res = " << solverPerfU.initialResidual()
	   << ", final rel res = " << relResU << nl
	   << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
	   << ", ClockTime = " << runTime.elapsedClockTime() << " s" 
	   << endl;
	
#       include "calculateEpsilonSigma.H"
#       include "writeFields.H"

        Info<< "ExecutionTime = "
            << runTime.elapsedCpuTime()
            << " s\n\n" << endl;
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
