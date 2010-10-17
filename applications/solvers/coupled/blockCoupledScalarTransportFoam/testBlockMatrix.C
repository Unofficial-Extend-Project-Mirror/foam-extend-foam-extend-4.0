// The FOAM Project // File: testCoeffField.C
/*
-------------------------------------------------------------------------------
 =========         | Application
 \\      /         |
  \\    /          | Name:   testCoeffField
   \\  /           | Family: Utility
    \\/            |
    F ield         | FOAM version: 2.2
    O peration     |
    A and          | Copyright (C) 1991-2003 Nabla Ltd.
    M anipulation  |          All Rights Reserved.
-------------------------------------------------------------------------------
DESCRIPTION
    Test coeff field and block matrix

AUTHOR
    Hrvoje Jasak

-------------------------------------------------------------------------------
*/

#include "argList.H"
#include "fieldTypes.H"
#include "blockLduMatrices.H"
#include "blockLduSolvers.H"
#include "Time.H"
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
            mesh.solver("HrvsVar")
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
