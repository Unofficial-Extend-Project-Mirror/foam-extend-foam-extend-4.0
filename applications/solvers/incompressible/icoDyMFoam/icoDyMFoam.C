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
    icoDyMFoam

Description
    Transient solver for incompressible, laminar flow of Newtonian fluids
    with dynamic mesh.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createDynamicFvMesh.H"
#   include "initContinuityErrs.H"
#   include "initTotalVolume.H"
#   include "createFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
#       include "readControls.H"
#       include "checkTotalVolume.H"

        // Make the fluxes absolute
        fvc::makeAbsolute(phi, U);

#       include "CourantNo.H"
#       include "setDeltaT.H"

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

        for (int corr = 0; corr < nCorr; corr++)
        {
            rAU = 1.0/UEqn.A();

            U = rAU*UEqn.H();
            phi = (fvc::interpolate(U) & mesh.Sf());
              //+ fvc::ddtPhiCorr(rAU, U, phi);

            adjustPhi(phi, U, p);

            for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
            {
                fvScalarMatrix pEqn
                (
                    fvm::laplacian(rAU, p) == fvc::div(phi)
                );

                pEqn.setReference(pRefCell, pRefValue);

                if (corr == nCorr - 1 && nonOrth == nNonOrthCorr)
                {
                    pEqn.solve(mesh.solutionDict().solver(p.name() + "Final"));
                }
                else
                {
                    pEqn.solve(mesh.solutionDict().solver(p.name()));
                }

                if (nonOrth == nNonOrthCorr)
                {
                    phi -= pEqn.flux();
                }
            }

#           include "continuityErrs.H"

            // Make the fluxes relative to the mesh motion
            fvc::makeRelative(phi, U);

            U -= rAU*fvc::grad(p);
            U.correctBoundaryConditions();
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
