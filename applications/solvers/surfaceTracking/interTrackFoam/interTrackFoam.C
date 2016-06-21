/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
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
    interfaceTrackinFoam

Description
    Incompressible laminar CFD code for simulation of a single bubble rising
    in a stil liquid. Interface between fluid phases is tracked using moving
    mesh.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "freeSurface.H"
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

    Info << "\nStarting time loop\n" << endl;

    for (runTime++; !runTime.end(); runTime++)
    {
        Info << "Time = " << runTime.value() << endl << endl;

        interface.moveMeshPointsForOldFreeSurfDisplacement();

        interface.updateDisplacementDirections();

        interface.predictPoints();

        Info<< "\nMax surface Courant Number = "
            << interface.maxCourantNumber() << endl << endl;

        while (pimple.loop())
        {
            // Update interface bc
            interface.updateBoundaryConditions();

            // Make the fluxes relative
            phi -= fvc::meshPhi(rho, U);

#           include "CourantNo.H"

            fvVectorMatrix UEqn
            (
                fvm::ddt(rho, U)
              + fvm::div(fvc::interpolate(rho)*phi, U, "div(phi,U)")
              - fvm::laplacian(mu, U)
            );

            solve(UEqn == -fvc::grad(p));

            // --- PISO loop
            while (pimple.correct())
            {
                volScalarField AU = UEqn.A();

                U = UEqn.H()/AU;

                phi = (fvc::interpolate(U) & mesh.Sf());

#               include "scalePhi.H"

                // Non-orthogonal pressure corrector loop
                while (pimple.correctNonOrthogonal())
                {
                    fvScalarMatrix pEqn
                    (
                        fvm::laplacian(1.0/AU, p)
                     == fvc::div(phi)
                    );

#                   include "setReference.H"

                    pEqn.solve();

                    if (pimple.finalNonOrthogonalIter())
                    {
                        phi -= pEqn.flux();
                    }
                }

#               include "continuityErrs.H"

                // Momentum corrector
                U -= fvc::grad(p)/AU;
                U.correctBoundaryConditions();
            }

            interface.correctPoints();

#           include "freeSurfaceContinuityErrs.H"
        }

#       include "volContinuity.H"

        Info << "Total surface tension force: "
            << interface.totalSurfaceTensionForce() << endl;

        vector totalForce =
            interface.totalViscousForce()
          + interface.totalPressureForce();

        Info << "Total force: " << totalForce << endl;

        runTime.write();

        Info << "ExecutionTime = "
            << scalar(runTime.elapsedCpuTime())
            << " s\n" << endl << endl;
    }

    Info << "End\n" << endl;

    return(0);
}

// ************************************************************************* //
