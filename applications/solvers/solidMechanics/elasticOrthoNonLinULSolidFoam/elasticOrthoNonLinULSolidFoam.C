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
    elasticOrthoGenDirULSolidFoam

Description
    Transient/steady-state segregated finite-volume solver for large strain
    elastic orthotropic solid bodies allowing for general principal material
    directions.

    Displacement increment field DU is solved for using an updated Lagrangian
    approach,
    also generating the Almansi strain tensor field epsilon and Cauchy stress
    tensor field sigma.

    At the end of each time-step, the mesh is moved and sigma, epsilon and C
    are rotated to the new configuration.

    Please cite:
    Cardiff P, Karac A & Ivankovic A, A Large Strain Finite Volume Method for
    Orthotropic Bodies with General Material Orientations, Computer Methods
    in Applied Mechanics & Engineering, Sep 2013,
    http://dx.doi.org/10.1016/j.cma.2013.09.008

Author
    Philip Cardiff UCD

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "constitutiveModel.H"
#include "transformField.H"
#include "transformGeometricField.H"
#include "pointPatchInterpolation.H"
#include "primitivePatchInterpolation.H"
#include "pointFields.H"
#include "twoDPointCorrector.H"
#include "leastSquaresVolPointInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#    include "setRootCase.H"
#    include "createTime.H"
#    include "createMesh.H"
#    include "createFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while(runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

#       include "readSolidMechanicsControls.H"

        int iCorr = 0;
        lduMatrix::solverPerformance solverPerf;
        scalar initialResidual = 1.0;
        lduMatrix::debug = 0;

        //- div(sigmaOld) should be zero but I will include
        //- it to make sure errors don't accumulate
        volVectorField* oldErrorPtr = NULL;
        if (ensureTotalEquilibrium)
        {
            oldErrorPtr = new volVectorField
            (
                fvc::d2dt2(rho.oldTime(), U.oldTime())
              - fvc::div(sigma)
            );
        }

        do
        {
            DU.storePrevIter();

            //- Updated lagrangian momentum equation
            fvVectorMatrix DUEqn
            (
                fvm::d2dt2(rho, DU)
              + fvc::d2dt2(rho, U)
              ==
                fvm::laplacian(K, DU, "laplacian(K,DU)")
              + fvc::div
                (
                    DSigma
                  - (K & gradDU)
                  + ( (sigma + DSigma) & gradDU ),
                    "div(sigma)"
                )
                //- fvc::laplacian(K, DU)
            );

            if (ensureTotalEquilibrium)
            {
                //- to stop accumulation of errors
                DUEqn += *oldErrorPtr;
            }

            solverPerf = DUEqn.solve();

            if (iCorr == 0)
            {
                initialResidual = solverPerf.initialResidual();
            }

            DU.relax();

            gradDU = fvc::grad(DU);

            //- for 2-D plane stress simulations, the zz component of gradDU
            //- ensures sigma.zz() is zero
            //- it is assumed that z is the empty direction
            //#         include "checkPlaneStress.H"

            //- sigma needs to be calculated inside the momentum loop as
            //- it is used in the momentum equation
            DEpsilon = symm(gradDU) + 0.5*symm(gradDU & gradDU.T());
            DSigma = C && DEpsilon;

            if (iCorr % infoFrequency == 0)
            {
                Info<< "\tTime " << runTime.value()
                    << ", Corr " << iCorr
                    << ", Solving for " << DU.name()
                    << " using " << solverPerf.solverName()
                    << ", res = " << solverPerf.initialResidual()
                    //<< ", rel res = " << relativeResidual
                    << ", inner iters " << solverPerf.nIterations() << endl;
            }
        }
        while
        (
            solverPerf.initialResidual() > convergenceTolerance
         && ++iCorr < nCorr
        );

        Info<< nl << "Time " << runTime.value() << ", Solving for " << DU.name()
            << ", Initial residual = " << initialResidual
            << ", Final residual = " << solverPerf.initialResidual()
            << ", No outer iterations " << iCorr
            << nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << endl;

#       include "moveMeshLeastSquares.H"
#       include "rotateFields.H"
#       include "writeFields.H"

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s\n\n"
            << endl;
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
