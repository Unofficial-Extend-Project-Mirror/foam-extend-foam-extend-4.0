/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
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
    sonicDyMFoam

Description
    Transient solver for trans-sonic/supersonic, laminar flow of a
    compressible gas with support for mesh motion and topological changes

    Updated from sonicFoamAutoMotion by Hrvoje Jasak

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "specie.H"
#include "basicPsiThermo.H"
#include "turbulenceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"

#   include "createTime.H"
#   include "createDynamicFvMesh.H"
#   include "createFields.H"
#   include "initContinuityErrs.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
#       include "readControls.H"
#       include "readFieldBounds.H"
#       include "compressibleCourantNo.H"
#       include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        bool meshChanged = mesh.update();

#       include "volContinuity.H"

        if (checkMeshCourantNo)
        {
#           include "meshCourantNo.H"
        }

        // Mesh motion update
//         if (correctPhi && meshChanged)
//         {
// #           include "correctPhi.H"
//         }

        if (meshChanged)
        {
#           include "CourantNo.H"
        }

        // --- PIMPLE loop
        label oCorr = 0;
        do
        {
            // Under-relax pDivU term
            pDivU.storePrevIter();
            pDivU = p*fvc::div(phi/fvc::interpolate(rho));
            pDivU.relax();

#           include "rhoEqn.H"
#           include "eEqn.H"
#           include "UEqn.H"

            // --- PISO loop
            for (int corr = 0; corr < nCorr; corr++)
            {
                U = UEqn.H()/UEqn.A();

#               include "limitU.H"

                surfaceScalarField phid
                (
                    "phid",
                    fvc::interpolate(psi)*
                    (
                        (fvc::interpolate(U) & mesh.Sf())
                      - fvc::meshPhi(rho, U)
                    )
                );

                for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
                {
                    // Store pressure for under-relaxation
                    p.storePrevIter();

                    fvScalarMatrix pEqn
                    (
                        fvm::ddt(psi, p)
                      + fvm::div(phid, p)
                      - fvm::laplacian(rho/UEqn.A(), p)
                    );

                    if
                    (
//                         oCorr == nOuterCorr - 1
                        corr == nCorr - 1
                     && nonOrth == nNonOrthCorr
                    )
                    {
                        pEqn.solve
                        (
                            mesh.solutionDict().solver(p.name() + "Final")
                        );
                    }
                    else
                    {
                        pEqn.solve(mesh.solutionDict().solver(p.name()));
                    }

                    if (nonOrth == nNonOrthCorr)
                    {
                        phi = pEqn.flux();
                    }


                    // Bound the pressure
                    if (min(p) < pMin || max(p) > pMax)
                    {
                        p.max(pMin);
                        p.min(pMax);
                        p.correctBoundaryConditions();
                    }

                    // Relax the pressure
                    p.relax();
                }

#               include "compressibleContinuityErrs.H"

                U -= fvc::grad(p)/UEqn.A();
                U.correctBoundaryConditions();

#               include "limitU.H"
            }

            // Recalculate density
            rho = thermo.rho();

            turbulence->correct();
        } while (++oCorr < nOuterCorr);

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
