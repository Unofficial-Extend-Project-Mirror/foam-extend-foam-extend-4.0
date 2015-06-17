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
    testBlockMatrix

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.

Description
    Test block matrix coefficient assembly

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "fieldTypes.H"
#include "blockLduMatrices.H"
#include "blockLduSolvers.H"
#include "foamTime.H"
#include "fvMesh.H"

#include "blockVector2Matrix.H"
#include "tensor2.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    blockScalarMatrix scalarMatrix(mesh);

    blockVectorMatrix blockMatrix(mesh);

    const label diagSize = mesh.lduAddr().size();
    const label ulSize = mesh.lduAddr().lowerAddr().size();
    const scalar diagCoeff = -2.0;
//     const scalar diagCoeff = -4.0;

//     Info << "Doing diagonal matrix" << endl;
    blockMatrix.upper() = scalarField(ulSize, 1.0);

    scalarField dScalar(diagSize, diagCoeff);
    dScalar[0] = -10000;
    dScalar[diagSize - 1] = -10000;

    blockMatrix.diag() = dScalar;

//     tensorField dTensor
//     (
//         diagSize,
//         tensor
//         (
//             diagCoeff,       0.0,  0.0,
//              0.0,      diagCoeff,  0.0,
//              0.0,            0.0, -1.0
//         )
//     );

//     dTensor[0] =
//         tensor
//         (
//             -10000.0,      0.0,  0.0,
//              0.0,     -10000.0,  0.0,
//              0.0,          0.0, -1.0
//         );

//     dTensor[diagSize - 1] =
//         tensor
//         (
//             -1.0,  0.0,  0.0,
//              0.0, -1.0,  0.0,
//              0.0,  0.0, -1.0
//         );

//     blockMatrix.diag() = dTensor;

    vectorField psi(diagSize, vector(0, 0, 0));
    vectorField source(diagSize, vector(0, 0, 0));
    source[0] = vector(0, 0, 0);
    source[diagSize - 1] = vector(10000, 0, 0);

//     psi[0] = vector(0, 0, 0);
//     psi[diagSize - 1] = vector(-1, 0, 0);

    BlockSolverPerformance<vector> solverPerf =
        blockVectorSolver::New
        (
            "HrvsVar",
            blockMatrix,
            mesh.solutionDict().solver("HrvsVar")
        )->solve(psi, source);

    Info << "Psi: " << psi << endl;
//     Info << "Psi: " << psi.component(vector::X) << endl;

    // Large block matrix
    BlockLduMatrix<vector2> vector2Matrix(mesh);

    vector2Matrix.diag().asScalar() =
        scalarField(vector2Matrix.diag().size(), 1);

    vector2Matrix.diag() +=
        Field<vector2>(vector2Matrix.diag().size(), vector2::one);

    vector2Matrix.upper() =
        Field<tensor2>(vector2Matrix.upper().size(), tensor2::one);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
