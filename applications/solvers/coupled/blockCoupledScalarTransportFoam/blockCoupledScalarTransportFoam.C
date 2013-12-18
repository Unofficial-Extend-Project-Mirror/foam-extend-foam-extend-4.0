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
    blockCoupledScalarTransportFoam

Description
    Solves two coupled transport equations in a block-coupled manner

        1) transport equation for a passive scalar
        2) diffusion only

    This resembles heat exchanging flow through a porous medium

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fieldTypes.H"
#include "Time.H"
#include "fvMesh.H"

#include "blockLduSolvers.H"
#include "VectorNFieldTypes.H"
#include "volVectorNFields.H"
#include "blockVectorNMatrices.H"

#include "blockMatrixTools.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
#   include "createFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating scalar transport\n" << endl;

#   include "CourantNo.H"

    for (runTime++; !runTime.end(); runTime++)
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

#       include "readSIMPLEControls.H"

        for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
        {
            fvScalarMatrix TEqn
            (
                fvm::div(phi, T)
              - fvm::laplacian(DT, T)
             ==
                 alpha*Ts
              - fvm::Sp(alpha, T)
            );

            TEqn.relax();

            fvScalarMatrix TsEqn
            (
              - fvm::laplacian(DTs, Ts)
             ==
                 alpha*T
              - fvm::Sp(alpha, Ts)
            );

            TsEqn.relax();

            // Prepare block system
            BlockLduMatrix<vector2> blockM(mesh);

            //- Transfer the coupled interface list for processor/cyclic/etc.
            // boundaries
            blockM.interfaces() = blockT.boundaryField().blockInterfaces();

            // Grab block diagonal and set it to zero
            Field<tensor2>& d = blockM.diag().asSquare();
            d = tensor2::zero;

            // Grab linear off-diagonal and set it to zero
            Field<vector2>& l = blockM.lower().asLinear();
            Field<vector2>& u = blockM.upper().asLinear();
            u = vector2::zero;
            l = vector2::zero;

            vector2Field& blockX = blockT.internalField();
            vector2Field blockB(mesh.nCells(), vector2::zero);

            //- Inset equations into block Matrix
            blockMatrixTools::insertEquation(0, TEqn, blockM, blockX, blockB);
            blockMatrixTools::insertEquation(1, TsEqn, blockM, blockX, blockB);

            //- Add off-diagonal terms and remove from block source
            forAll(d, i)
            {
                d[i](0, 1) = -alpha.value()*mesh.V()[i];
                d[i](1, 0) = -alpha.value()*mesh.V()[i];

                blockB[i][0] -= alpha.value()*blockX[i][1]*mesh.V()[i];
                blockB[i][1] -= alpha.value()*blockX[i][0]*mesh.V()[i];
            }

            //- Block coupled solver call
            BlockSolverPerformance<vector2> solverPerf =
                BlockLduSolver<vector2>::New
                (
                    blockT.name(),
                    blockM,
                    mesh.solutionDict().solver(blockT.name())
                )->solve(blockX, blockB);

            solverPerf.print();

            // Retrieve solution
            blockMatrixTools::blockRetrieve(0, T.internalField(), blockX);
            blockMatrixTools::blockRetrieve(1, Ts.internalField(), blockX);

            T.correctBoundaryConditions();
            Ts.correctBoundaryConditions();
        }

        runTime.write();
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
