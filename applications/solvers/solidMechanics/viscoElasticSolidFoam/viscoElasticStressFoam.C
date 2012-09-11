/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-6 H. Jasak All rights reserved
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
    viscoElasticStressedFoam

Description
    Transient/steady-state segregated finite-volume solver for small strain
    visco elastic solid bodies.

    Displacement increment field DU is solved for using a total Lagrangian
    approach, also generating the strain tensor field epsilon and stress
    tensor field sigma.


\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "rheologyModel.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
# include "setRootCase.H"

# include "createTime.H"

# include "createMesh.H"

# include "createFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

  Info<< "\nCalculating displacement field\n" << endl;

  lduMatrix::debug = 0;
  
  scalar m = 0.5;
  
  for (runTime++; !runTime.end(); runTime++)
    {
      Info<< "Time: " << runTime.timeName() << nl << endl;
      
#     include "readStressedFoamControls.H"

      volScalarField mu = 
	rheology.mu(m*runTime.deltaT().value());
      volScalarField lambda = 
	rheology.lambda(m*runTime.deltaT().value());
      
      Info << "mu = " << average(mu.internalField()) << endl;
      Info << "lambda = " << average(lambda.internalField()) << endl;
      
      int iCorr = 0;
      lduMatrix::solverPerformance solverPerf;
      scalar initialResidual = 0;
      scalar residual = GREAT;
      
      do
        {
	  DU.storePrevIter();
	  
	  fvVectorMatrix DUEqn
            (
	     fvm::d2dt2(rho,DU)
             ==
	     fvm::laplacian(2*mu+lambda, DU, "laplacian(DDU,DU)")
	     + fvc::div
	     (
	      mu*gradDU.T() 
	      + lambda*(I*tr(gradDU))
	      - (mu + lambda)*gradDU
	      + DSigmaCorr,
	      "div(sigma)"
	      )
	     );

	  solverPerf = DUEqn.solve();
	  
	  DU.relax();
	  
	  if(iCorr == 0)
            {
	      initialResidual = solverPerf.initialResidual();
            }
	  
	  gradDU = fvc::grad(DU);
	  
#         include "calculateDEpsilonDSigma.H"
            
        }
      while
        (
	 solverPerf.initialResidual() > convergenceTolerance
         && ++iCorr < nCorr
	 );
      
      Info << "Solving for " << DU.name() << " using "
	   << solverPerf.solverName() << " solver"
	   << ", Initial residual = " << initialResidual
	   << ", Final residual = " << solverPerf.initialResidual()
	   << ", No outer iterations " << iCorr
	   << ", Relative error: " << residual << endl;
      
      U += DU;

      epsilon += DEpsilon;
      
#     include "calculateSigmaDSigmaCorr.H"

#     include "writeFields.H"

      Info<< "ExecutionTime = "
	  << runTime.elapsedCpuTime()
	  << " s\n\n" << endl;
    }

  Info<< "End\n" << endl;
  
  return(0);
}


// ************************************************************************* //
