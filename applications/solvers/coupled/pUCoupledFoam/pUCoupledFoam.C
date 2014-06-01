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
    pUCoupledFoam

Description
    Steady-state solver for incompressible, turbulent flow, with implicit
    coupling between pressure and velocity achieved by BlockLduMatrix
    Turbulence is in this version solved using the existing turbulence
    structure.

Authors
    Klas Jareteg, Chalmers University of Technology,
    Vuko Vukcevic, FMENA Zagreb.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "RASModel.H"

#include "VectorNFieldTypes.H"
#include "volVectorNFields.H"
#include "blockLduSolvers.H"
#include "blockVectorNMatrices.H"

#include "blockMatrixTools.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
#   include "createFields.H"
#   include "initContinuityErrs.H"
#   include "readBlockSolverControls.H"

    // Calculate coupling matrices only once since the mesh doesn't change and
    // implicit div and grad operators are only dependant on Sf. Actually
    // coupling terms (div(U) and grad(p)) in blockMatrix do not change, so they
    // could be inserted only once, resetting other parts of blockMatrix to zero
    // at the end of each time step. VV, 30/April/2014
#   include "calculateCouplingMatrices.H"

    Info<< "\nStarting time loop\n" << endl;
    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        p.storePrevIter();

        // Initialize block matrix
#       include "initializeBlockMatrix.H"

        // Assemble and insert momentum equation
#       include "UEqn.H"

        // Assemble and insert pressure equation
#       include "pEqn.H"

        // Insert coupling, updating the boundary contributions
        // Last argument in insertBlockCoupling says if the first location
        // should be incremented. This is needed for arbitrary positioning
        // of U and p in the system. This could be better. VV, 30/April/2014
        blockMatrixTools::insertBlockCoupling(3, 0, UInp, U, A, b, false);
        blockMatrixTools::insertBlockCoupling(0, 3, pInU, p, A, b, true);

        // Solve the block matrix
        BlockSolverPerformance<vector4> solverPerf =
            BlockLduSolver<vector4>::New
            (
                word("Up"),
                A,
                mesh.solutionDict().solver("Up")
            )->solve(Up, b);

        solverPerf.print();

        // Retrieve solution
        blockMatrixTools::retrieveSolution(0, U.internalField(), Up);
        blockMatrixTools::retrieveSolution(3, p.internalField(), Up);

        U.correctBoundaryConditions();
        p.correctBoundaryConditions();

        phi = (fvc::interpolate(U) & mesh.Sf()) + pEqn.flux() + presSource;

#       include "continuityErrs.H"

        p.relax();

        turbulence->correct();
        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}
