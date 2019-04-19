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
    transientDyMFoam

Description
    Transient solver for incompressible, turbulent flow, with implicit
    coupling between pressure and velocity achieved by fvBlockMatrix.

    Turbulence is solved using the existing turbulence model structure.

    The solver supports dynamic mesh changes

Authors
    Hrvoje Jasak, Wikki Ltd.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvBlockMatrix.H"
#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"
#include "dynamicFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createDynamicFvMesh.H"
#   include "createFields.H"
#   include "initContinuityErrs.H"
#   include "initConvergenceCheck.H"
#   include "createControls.H"

    Info<< "\nStarting time loop\n" << endl;
    while (runTime.run())
    {
#       include "readBlockSolverControls.H"
#       include "readFieldBounds.H"
#       include "CourantNo.H"
#       include "setDeltaT.H"

        // Make the fluxes absolute
        fvc::makeAbsolute(phi, U);

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        bool meshChanged = mesh.update();
        reduce(meshChanged, orOp<bool>());

#       include "volContinuity.H"

        if (correctPhi && meshChanged)
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

        if (meshChanged)
        {
#           include "CourantNo.H"
        }

        for (label i = 0; i < nOuterCorrectors; i++)
        {
            p.storePrevIter();

            // Initialize the Up block system
            fvBlockMatrix<vector4> UpEqn(Up);

            // Assemble and insert momentum equation
#           include "UEqn.H"

            // Assemble and insert pressure equation
#           include "pEqn.H"

            // Assemble and insert coupling terms
#           include "couplingTerms.H"

            // Solve the block matrix
            residual = UpEqn.solve();
            maxResidual = cmptMax(residual.initialResidual());

            // Retrieve solution
            UpEqn.retrieveSolution(0, U.internalField());
            UpEqn.retrieveSolution(3, p.internalField());

            U.correctBoundaryConditions();
            p.correctBoundaryConditions();

            phi = (fvc::interpolate(U) & mesh.Sf())
                + pEqn.flux()
                + presSource;

            // Make the fluxes relative to the mesh motion
            fvc::makeRelative(phi, U);

#           include "movingMeshContinuityErrs.H"

#           include "boundPU.H"

            p.relax();

            turbulence->correct();
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

#       include "convergenceCheck.H"
    }

    Info<< "End\n" << endl;

    return 0;
}
