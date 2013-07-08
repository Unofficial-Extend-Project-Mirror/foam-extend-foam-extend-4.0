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


\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "rheologyModel.H"
#include "thermalModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
# include "setRootCase.H"

# include "createTime.H"

# include "createMesh.H"

# include "createFields.H"

# include "readSigmaExpMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

  Info<< "\nCalculating displacement field\n" << endl;

  while(runTime.loop())
    {
      Info<< "Time: " << runTime.timeName() << nl << endl;
      
#     include "readStressedFoamControls.H"
      
      int iCorr = 0;
      scalar initialResidual = GREAT;
      scalar residual = GREAT;
      lduMatrix::solverPerformance solverPerfU;
      lduMatrix::solverPerformance solverPerfT;
     
      lduMatrix::debug=0;

      do
        {
	  U.storePrevIter();
	 
#         include "calculateSigmaExp.H"

	  //- energy equation
	  fvScalarMatrix TEqn
	    (
	     fvm::ddt(rhoC, T) == fvm::laplacian(k, T, "laplacian(k,T)")
	     );

	  solverPerfT = TEqn.solve();

	  T.relax();

	  Info << "\tTime " << runTime.value()
	       << ", Corrector " << iCorr << nl
	       << "\t\tSolving for " << T.name()
	       << " using " << solverPerfT.solverName()
	       << ", residual = " << solverPerfT.initialResidual() << endl;

	  //- linear momentum equaiton
	  fvVectorMatrix UEqn
            (
	     fvm::d2dt2(rho, U)
	     ==
	     fvm::laplacian(2*mu + lambda, U, "laplacian(DU,U)")
	     + sigmaExp
	     - fvc::grad(threeKalpha*(T-T0),"grad(threeKalphaDeltaT)")
	     );

	  solverPerfU = UEqn.solve();

	  if(iCorr == 0)
	    {
	      initialResidual = max
		(
		 solverPerfU.initialResidual(),
		 solverPerfT.initialResidual()
		 );
	    }

	  residual = max
	    (
	     solverPerfU.initialResidual(),
	     solverPerfT.initialResidual()
	     );

	  U.relax();
	  
	  gradU = fvc::grad(U);
	  
	  Info << "\t\tSolving for " << U.name()
	       << " using " << solverPerfU.solverName()
	       << ", residual = " << solverPerfU.initialResidual() << endl;
        }
	while
	  (
	   residual > convergenceTolerance
	   &&
	   ++iCorr < nCorr
	   );
	
        Info << nl << "Time " << runTime.value()
	     << ", Solving for " << U.name() 
	     << ", Solving for " << T.name() 
	     << ", Initial residual = " << initialResidual 
	     << ", Final U residual = " << solverPerfU.initialResidual()
	     << ", Final T residual = " << solverPerfT.initialResidual()
	     << ", No outer iterations " << iCorr
	     << nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
	     << "  ClockTime = " << runTime.elapsedClockTime() << " s" 
	     << endl;
	
	lduMatrix::debug=0;

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
