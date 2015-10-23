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
    MRFPorousFoam

Description
    Steady-state solver for incompressible, turbulent flow, with implicit
    coupling between pressure and velocity achieved by fvBlockMatrix.
    Turbulence is in this version solved using the existing turbulence
    structure.

    Added support for MRF and porous zones

Authors
    Hrvoje Jasak, Wikki Ltd.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvBlockMatrix.H"
#include "singlePhaseTransportModel.H"
#include "RASModel.H"
#include "MRFZones.H"
#include "porousZones.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
#   include "createFields.H"
#   include "createMRF.H"
#   include "createPorous.H"
#   include "initContinuityErrs.H"
#   include "initConvergenceCheck.H"

    Info<< "\nStarting time loop\n" << endl;
    while (runTime.loop())
    {
#       include "readBlockSolverControls.H"
#       include "readFieldBounds.H"

        Info<< "Time = " << runTime.timeName() << nl << endl;

        p.storePrevIter();

        // Initialize the Up block system (matrix, source and reference to Up)
        fvBlockMatrix<vector4> UpEqn(Up);

        // Assemble and insert momentum equation
#       include "UEqn.H"

        // Assemble and insert pressure equation
#       include "pEqn.H"

        // Assemble and insert coupling terms
#       include "couplingTerms.H"

        // Solve the block matrix
        maxResidual = cmptMax(UpEqn.solve().initialResidual());

        // Retrieve solution
        UpEqn.retrieveSolution(0, U.internalField());
        UpEqn.retrieveSolution(3, p.internalField());

        U.correctBoundaryConditions();
        p.correctBoundaryConditions();

        phi = (fvc::interpolate(U) & mesh.Sf()) + pEqn.flux() + presSource;

        // Make flux relative in rotating zones
        mrfZones.relativeFlux(phi);

#       include "continuityErrs.H"

#       include "boundPU.H"

        p.relax();

        turbulence->correct();
        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

#       include "convergenceCheck.H"
    }

    Info<< "End\n" << endl;

    return 0;
}
