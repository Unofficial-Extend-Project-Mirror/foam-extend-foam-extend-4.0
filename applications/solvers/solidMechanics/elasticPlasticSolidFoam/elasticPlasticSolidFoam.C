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
    elasticPlasticSolidFoam

Description
    Transient/steady-state segregated finite-volume solver for small strain
    elastic plastic solid bodies.

    Displacement increment field DU is solved for using a total Lagrangian
    approach, also generating the strain tensor field epsilon, the plastic
    strain field epsilonP and stress tensor field sigma.

Author
    A. Karac, A. Ivankovic,
    P. Cardiff UCD

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
//#include "plasticityModel.H"
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
  
  Info<< "\nStarting time loop\n" << endl;
  
  for (runTime++; !runTime.end(); runTime++)
    {
      Info<< "Time: " << runTime.timeName() << nl << endl;
      
#     include "readStressedFoamControls.H"
      
      int iCorr = 0;
      scalar initialResidual = 0;
      scalar relativeResidual = 1.0;
      scalar plasticResidual = 1.0;
      lduMatrix::solverPerformance solverPerf;
      lduMatrix::debug = 0;
      
      const volSymmTensorField& DEpsilonP = rheology.DEpsilonP();
      
//       volVectorField* oldErrorPtr = NULL;
//       if(ensureTotalEquilibrium)
// 	{
	  //const volScalarField& beta =
	  //mesh.objectRegistry::lookupObject<volScalarField>("beta");
// 	  oldErrorPtr = new volVectorField
// 	    (
// 	     rho*fvc::d2dt2(U.oldTime())
// // 	     - fvc::div(sigma.oldTime())
// 	     - fvc::div(mesh.Sf() & fvc::interpolate(sigma.oldTime()))
// 	     );
// 	}

      do
	{
	  DU.storePrevIter();
	  DEpsilonP.storePrevIter();

#         include "calculateDivDSigmaExp.H"
	  	  
	  fvVectorMatrix DUEqn
	    (
	     rho*fvm::d2dt2(DU)
	     ==
	     fvm::laplacian(2*muf + lambdaf, DU, "laplacian(DDU,DU)")
	     + divDSigmaExp
	     - fvc::div(2*muf*(mesh.Sf() & fvc::interpolate(DEpsilonP)))
	     );

	  // if(ensureTotalEquilibrium)
	  //   {
	      // 	      const volScalarField& beta =
	      // 		mesh.objectRegistry::lookupObject<volScalarField>("beta");
	      // 	      oldErrorPtr = new volVectorField
	      // 		(
	      // 		 rho*fvc::d2dt2(U.oldTime())
	      // 		 - fvc::div(sigma.oldTime())
	      // 		 );
	      // 	      DUEqn += *oldErrorPtr;
	      //DUEqn -= fvc::div(mesh.Sf() & fvc::interpolate(sigma, "sigma"));
	    // }

	  solverPerf = DUEqn.solve();
	  
	  if(iCorr == 0)
            {
	      initialResidual = solverPerf.initialResidual();
            }
	  
	  if(aitkenRelax)
	    {
#             include "aitkenRelaxation.H"
	    }
	  else
	    {	  
	      DU.relax();
	    }
	  
	  gradDU = fvc::grad(DU);

#         include "calculateRelativeResidual.H"
#         include "calculateDEpsilonDSigma.H"

	  // correct plasticity
	  rheology.correct();

	  // update mu and lambda for non-linear elastic
	  //mu = rheology.newMu();
	  //lambda = rheology.newLambda();
	  //muf = fvc::interpolate(mu);
	  //lambdaf = fvc::interpolate(lambda);

          if(iCorr % infoFrequency == 0)
            {
              Info << "\tTime " << runTime.value()
                   << ", Corrector " << iCorr
                   << ", Solving for " << DU.name()
		   << " using " << solverPerf.solverName()
                   << ", res = " << solverPerf.initialResidual()
                   << ", rel res = " << relativeResidual
                   << ", plastic res = " << plasticResidual;
	      if(aitkenRelax) Info << ", aitken = " << aitkenTheta;
	      Info << ", inner iters = " << solverPerf.nIterations() << endl;
            }
	}
      while
	(
	 iCorr++ < 2
	 ||
	 (//solverPerf.initialResidual() > convergenceTolerance
	  relativeResidual > convergenceTolerance
	  &&
	  iCorr < nCorr)
	 );
      
      Info << nl << "Time " << runTime.value() << ", Solving for " << DU.name() 
	   << ", Initial residual = " << initialResidual 
	   << ", Final residual = " << solverPerf.initialResidual()
	   << ", Final rel residual = " << relativeResidual
	   << ", No outer iterations " << iCorr << endl;
      
      lduMatrix::debug = 1;
      
      // update total quantities
      U += DU;
      epsilon += DEpsilon;
      epsilonP += rheology.DEpsilonP();
      sigma += DSigma;

      // update yields stresses
      rheology.updateYieldStress();
      
#     include "writeFields.H"
#     include "writeHistory.H"
      
      Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
	  << "  ClockTime = " << runTime.elapsedClockTime() << " s"
	  << nl << endl;
    }
  
  Info<< "End\n" << endl;
  
  return(0);
}


// ************************************************************************* //
