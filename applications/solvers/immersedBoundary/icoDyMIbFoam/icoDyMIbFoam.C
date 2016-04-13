/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
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
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    icoDyMOversetFoam

Description
    Transient solver for incompressible, laminar flow of Newtonian fluids
    with dynamic mesh and immersed boundary mesh support.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "immersedBoundaryFvPatch.H"
#include "immersedBoundaryAdjustPhi.H"
#include "pimpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createDynamicFvMesh.H"

    pimpleControl pimple(mesh);

#   include "createFields.H"
#   include "initContinuityErrs.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Make the fluxes absolute
        fvc::makeAbsolute(phi, U);

        bool meshChanged = mesh.update();
        reduce(meshChanged, orOp<bool>());
        Info<< "Mesh update" << meshChanged << endl;
#       include "createIbMasks.H"

        // Make the fluxes relative to the mesh motion
        fvc::makeRelative(phi, U);

#       include "CourantNo.H"

        // Pressure-velocity corrector
        while (pimple.loop())
        {
            fvVectorMatrix UEqn
            (
                fvm::ddt(U)
              + fvm::div(phi, U)
              - fvm::laplacian(nu, U)
            );

            if (pimple.momentumPredictor())
            {
                solve(UEqn == -fvc::grad(p));
            }

            // --- PISO loop
            while (pimple.correct())
            {
                volScalarField rUA = 1.0/UEqn.A();

                U = rUA*UEqn.H();
                // Immersed boundary update
                U.correctBoundaryConditions();

                phi = faceIbMask*(fvc::interpolate(U) & mesh.Sf());

                // Adjust immersed boundary fluxes
                immersedBoundaryAdjustPhi(phi, U);
                adjustPhi(phi, U, p);

                // Non-orthogonal pressure corrector loop
                while (pimple.correctNonOrthogonal())
                {
                    fvScalarMatrix pEqn
                    (
                        fvm::laplacian(rUA, p) == fvc::div(phi)
                    );

                    pEqn.setReference(pRefCell, pRefValue);
                    pEqn.solve
                    (
                        mesh.solutionDict().solver
                        (
                            p.select(pimple.finalInnerIter())
                        )
                    );

                    if (pimple.finalNonOrthogonalIter())
                    {
                        phi -= pEqn.flux();
                    }
                }

#               include "immersedBoundaryContinuityErrs.H"

                // Make the fluxes relative to the mesh motion
                fvc::makeRelative(phi, U);

                U -= rUA*fvc::grad(p);
                U.correctBoundaryConditions();
            }
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
