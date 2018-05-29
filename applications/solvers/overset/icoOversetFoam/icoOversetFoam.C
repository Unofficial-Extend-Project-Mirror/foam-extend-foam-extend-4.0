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
    icoOversetFoam

Description
    Transient solver for incompressible, laminar flow of Newtonian fluids
    with overset mesh support.

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pisoControl.H"

#include "oversetMesh.H"
#include "oversetFvPatchFields.H"
#include "oversetAdjustPhi.H"
#include "globalOversetAdjustPhi.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    pisoControl piso(mesh);

#   include "createOversetMasks.H"
#   include "createFields.H"
#   include "initContinuityErrs.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

#       include "oversetCourantNo.H"

        // Momentum equation
        fvVectorMatrix UEqn
        (
            fvm::ddt(U)
          + fvm::div(phi, U)
          - fvm::laplacian(nu, U)
        );

        if (piso.momentumPredictor())
        {
            solve(UEqn == -fvc::grad(p));
        }

        // --- PISO loop

        while (piso.correct())
        {
            p.boundaryField().updateCoeffs();

            rAU = 1.0/UEqn.A();
            oversetFvPatchScalarField::oversetInterpolate(rAU); // Overset update

            U = rAU*UEqn.H();
            oversetFvPatchVectorField::oversetInterpolate(U); // Overset update

            phi = fvc::interpolate(U) & mesh.Sf();

            // Adjust overset fluxes
            oversetAdjustPhi(phi, U); // Fringe flux adjustment
            globalOversetAdjustPhi(phi, U, p); // Global flux adjustment

            while (piso.correctNonOrthogonal())
            {
                fvScalarMatrix pEqn
                (
                    fvm::laplacian(rAU, p) == fvc::div(phi)
                );

                // Adjust non-orthogonal fringe fluxes if necessary
                om.correctNonOrthoFluxes(pEqn, U);

                pEqn.setReference(pRefCell, pRefValue);
                pEqn.solve
                (
                    mesh.solutionDict().solver(p.select(piso.finalInnerIter()))
                );

                if (piso.finalNonOrthogonalIter())
                {
                    phi -= pEqn.flux();
                }

                // Perform overset interpolation (after flux reconstruction)
                oversetFvPatchScalarField::oversetInterpolate(p);
            }

#           include "oversetContinuityErrs.H"

            U -= rAU*fvc::grad(p);
            U.correctBoundaryConditions();
            oversetFvPatchVectorField::oversetInterpolate(U); // Overset update
            // Note: if implicit fringe is switched on, we are doing the
            // interpolation twice (once in correctBoundaryConditions and once
            // in oversetInterpolate). Reorganize. VV, 4/Oct/2016.

            // Make the fluxes relative to the mesh motion
            fvc::makeRelative(phi, U);
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
