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
    Transient/steady-state segregated finite-volume solver for large strain
    elastic plastic solid bodies.

    Displacement increment field DU is solved for using an updated Lagrangian
    approach, also generating the strain field epsilon, the plastic strain
    field epsilonP and the stress tensor field sigma.

    With optional multi-material solid interface correction ensuring
    correct tractions on multi-material interfaces, HOWEVER, tractions
    on interface will be incorrect when there is plasticity or large strain
    in the interface cells. Correction needs to be derived for plasticity
    and large strain.

Author
    Philip Cardiff
    multi-material by Tukovic et al. 2012

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "plasticityModel.H"
#include "solidInterface.H"
#include "volPointInterpolation.H"
#include "pointPatchInterpolation.H"
#include "primitivePatchInterpolation.H"
#include "twoDPointCorrector.H"
#include "pointFields.H"
#include "plane.H"
#include "meshSearch.H"
#include "leastSquaresVolPointInterpolation.H"

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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info << "\nStarting time loop\n" << endl;

    for (runTime++; !runTime.end(); runTime++)
    {
        Info<< "Time: " << runTime.timeName() << nl << endl;

#       include "readStressedFoamControls.H"

        int iCorr = 0;
        lduMatrix::solverPerformance solverPerf;
        scalar initialResidual = 0;
        scalar relativeResidual = GREAT;
        lduMatrix::debug = 0;

        const volSymmTensorField& DEpsilonP = rheology.DEpsilonP();

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
                fvm::d2dt2(rho, DU)
              ==
                fvm::laplacian(2*muf + lambdaf, DU, "laplacian(DDU,DU)")
              + divDSigmaExp
              + divDSigmaLargeStrainExp
              - fvc::div(2*muf*(mesh.Sf() & fvc::interpolate(DEpsilonP)))
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

            DF = gradDU.T();

#           include "calculateRelativeResidual.H"

            rheology.correct();
            mu = rheology.newMu();
            lambda = rheology.newLambda();
            muf = fvc::interpolate(rheology.newMu());
            lambdaf = fvc::interpolate(rheology.newLambda());
            if(solidInterfaceCorr)
            {
                solidInterfacePtr->modifyProperties(muf, lambdaf);
            }

#           include "calculateDEpsilonDSigma.H"

            Info << "\tTime " << runTime.value()
                << ", Corrector " << iCorr
                << ", Solving for " << DU.name()
                << " using " << solverPerf.solverName()
                << ", residual = " << solverPerf.initialResidual()
                << ", relative residual = " << relativeResidual << endl;
        }
        while
        (
            //relativeResidual
            solverPerf.initialResidual() > convergenceTolerance
            && ++iCorr < nCorr
        );

        Info << nl << "Time " << runTime.value() << ", Solving for " << DU.name()
            << ", Initial residual = " << initialResidual
            << ", Final residual = " << solverPerf.initialResidual()
            << ", No outer iterations " << iCorr << endl;

        lduMatrix::debug = 1;

        U += DU;

        epsilon += DEpsilon;

        epsilonP += DEpsilonP;

        volSymmTensorField DEpsilonE = DEpsilon - DEpsilonP;

        epsilonE += DEpsilonE;

        sigma += DSigma;

        rheology.updateYieldStress();

#       include "rotateFields.H"

#       include "moveMesh.H"

#       include "writeFields.H"

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
