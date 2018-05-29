/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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
    icoDyMFoam

Description
    Transient solver for incompressible, laminar flow of Newtonian fluids
    with dynamic mesh.
    Consistent formulation without time-step and relaxation dependence by Jasak
    and Tukovic.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "pisoControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createDynamicFvMesh.H"

    pisoControl piso(mesh);

#   include "initContinuityErrs.H"
#   include "initTotalVolume.H"
#   include "createFields.H"
#   include "createControls.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
#       include "readControls.H"
#       include "checkTotalVolume.H"
#       include "CourantNo.H"
#       include "setDeltaT.H"

        // Make the fluxes absolute
        fvc::makeAbsolute(phi, U);

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        bool meshChanged = mesh.update();
        reduce(meshChanged, orOp<bool>());

#       include "volContinuity.H"

        if (correctPhi && (mesh.moving() || meshChanged))
        {
            // Fluxes will be corrected to absolute velocity
            // HJ, 6/Feb/2009
#           include "correctPhi.H"
        }

        // Make the fluxes relative to the mesh motion
        fvc::makeRelative(phi, U);

        if (mesh.moving() && checkMeshCourantNo)
        {
#           include "meshCourantNo.H"
        }

#       include "UEqn.H"

        // --- PISO loop

        // Prepare clean 1/a_p without time derivative contribution
        rAU = 1.0/HUEqn.A();

        while (piso.correct())
        {
            // Calculate U from convection-diffusion matrix
            U = rAU*HUEqn.H();

            // Consistently calculate flux
            piso.calcTransientConsistentFlux(phi, U, rAU, ddtUEqn);

            adjustPhi(phi, U, p);

            // Non-orthogonal pressure corrector loop
            while (piso.correctNonOrthogonal())
            {
                fvScalarMatrix pEqn
                (
                    fvm::laplacian
                    (
                        fvc::interpolate(rAU)/piso.aCoeff(U.name()),
                        p,
                        "laplacian(rAU," + p.name() + ')'
                    )
                 ==
                    fvc::div(phi)
                );

                pEqn.setReference(pRefCell, pRefValue);
                pEqn.solve
                (
                    mesh.solutionDict().solver(p.select(piso.finalInnerIter()))
                );

                if (piso.finalNonOrthogonalIter())
                {
                    phi -= pEqn.flux();
                }
            }

#           include "movingMeshContinuityErrs.H"

            // Consistently reconstruct velocity after pressure equation.
            // Note: flux is made relative inside the function
            piso.reconstructTransientVelocity(U, phi, ddtUEqn, rAU, p);
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
