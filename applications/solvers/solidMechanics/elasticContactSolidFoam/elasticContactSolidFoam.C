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
    elasticContactSolidFoam

Description
    Transient/steady-state segregated finite-volume solver for small strain
    elastic solid bodies in contact, using an total strain total Lagrangian
    approach.

    Works in parallel but mesh.movePoints sometimes fails for some unknown
    reason depending on the decomposition.

    Solves for the displacement increment vector field DU, also generating the
    stress tensor field sigma.

    It is only for frictionless contact, friction not implemented yet.

Author
    Philip Cardiff

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "rheologyModel.H"
#include "contactProblem.H"

#include "volPointInterpolation.H"
#include "pointPatchInterpolation.H"
#include "primitivePatchInterpolation.H"
#include "fixedValuePointPatchFields.H"
#include "pointFields.H"
#include "pointMesh.H"
#include "pointBoundaryMesh.H"
#include "primitivePatchInterpolation.H"
#include "twoDPointCorrector.H"
#include "plane.H"
#include "meshSearch.H"
//#include "leastSquaresVolPointInterpolation.H"

#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
# include "setRootCase.H"

# include "createTime.H"

# include "createMesh.H"

# include "createFields.H"

# include "readDivSigmaExpMethod.H"

# include "createGlobalToLocalFaceZonePointMap.H"

  // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
  
  Info<< "\nStarting time loop\n" << endl;
  
  for (runTime++; !runTime.end(); runTime++)
    {
      Info<< "Time: " << runTime.timeName() << endl;
      
#     include "readContactControls.H"

#     include "readStressedFoamControls.H"
  
      //-- for moving the mesh and then back again
      vectorField oldMeshPoints = mesh.allPoints();
      
      int iCorr = 0;
      lduMatrix::solverPerformance solverPerf;
      word solverName;
      lduMatrix::debug = 0;
      scalar residual = GREAT;
      scalar initialResidual = 0;
      scalar relativeResidual = GREAT;

      //- Predictor step
      if (predictor)
        {
	  Info << "\nPredicting U, gradU and snGradU based on V, gradV and snGradV\n" << endl;
	  U += V*runTime.deltaT();
	  gradU += gradV*runTime.deltaT();
	  snGradU += snGradV*runTime.deltaT();
        }

      do //- start of momentum loop
	{
	  U.storePrevIter();
	  
	  //- correct the contact boundaries 
          if(iCorr % uEqnContactCorrFreq == 0)
	    {
	      Info << "\t\tCorrecting contact in the momentum loop "
		   << "iteration: " << iCorr
		   << ", residual: " << residual
		   << endl;
	      //#                 include "moveMeshLeastSquares.H"
#             include "moveSolidMesh.H"
	      contact.correct();
	      mesh.movePoints(oldMeshPoints);
	    }
  
#         include "calculateDivSigmaExp.H"

	  fvVectorMatrix UEqn
	    (
	     fvm::d2dt2(rho, U)
	     ==
	     fvm::laplacian(2*mu + lambda, U, "laplacian(DU,U)")
	     + divSigmaExp
	     );
	  
	  solverPerf = UEqn.solve();
	  
	  U.relax();
	  
	  solverName = solverPerf.solverName();
	  
	  gradU = fvc::grad(U);
	  snGradU = fvc::snGrad(U);

	  residual = solverPerf.initialResidual();
	  
	  if(iCorr == 0)
	    {
	      initialResidual = solverPerf.initialResidual();
	    }
	  
#         include "calculateRelativeResidual.H"	      

	  Info << "\tTime " << runTime.value()
	       << ", Corrector " << iCorr
	       << ", Solving for " << U.name()
	       << " using " << solverPerf.solverName()
	       << ", residual = " << solverPerf.initialResidual()
	       << ", relative residual = " << relativeResidual << endl;
	} //- end of momentum loop
      while
	(
	 //relativeResidual > convergenceTolerance
	 residual > convergenceTolerance
	 &&
	 ++iCorr < nCorr
	 );
      
      // Print out info per contact iteration
      Info << "\t\tSolving for " << U.name()
	   << " using " << solverName
	   << ", Initial residual = " << initialResidual
	   << ", Final residual = " << solverPerf.initialResidual()
	   << ", No outer iterations " << iCorr << endl;
      
      lduMatrix::debug = 1;
      
      V = fvc::ddt(U);
      gradV = fvc::ddt(gradU);
      snGradV = (snGradU - snGradU.oldTime())/runTime.deltaT();

#     include "calculateEpsilonSigma.H"
      
#     include "writeFields.H"
        
      //#     include "moveMeshLeastSquares.H"
      //#     include "moveSolidMesh.H"
      //#     include "printContactResults.H"
      //mesh.movePoints(oldMeshPoints);
      
      Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
	  << "  ClockTime = " << runTime.elapsedClockTime() << " s"
	  << endl << endl;
    }
  
  Info<< "End\n" << endl;
  
  return(0);
}


// ************************************************************************* //
