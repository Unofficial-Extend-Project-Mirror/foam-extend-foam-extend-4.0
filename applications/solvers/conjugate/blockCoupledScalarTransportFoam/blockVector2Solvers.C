/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-6 H. Jasak All rights reserved
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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

\*---------------------------------------------------------------------------*/

#include "vector2Field.H"
#include "tensor2Field.H"
#include "ExpandTensorN.H"
#include "ExpandTensorNField.H"
#include "blockLduMatrices.H"
#include "addToRunTimeSelectionTable.H"

#include "blockLduPrecons.H"
#include "BlockNoPrecon.H"
#include "blockDiagonalPrecons.H"
#include "blockGaussSeidelPrecons.H"
#include "BlockCholeskyPrecon.H"

#include "blockLduSmoothers.H"
#include "blockGaussSeidelSmoothers.H"
#include "BlockILUSmoother.H"

#include "blockLduSolvers.H"
#include "BlockDiagonalSolver.H"
#include "BlockBiCGStabSolver.H"
#include "BlockCGSolver.H"
#include "BlockGaussSeidelSolver.H"
#include "BlockGMRESSolver.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// Preconditioners
typedef BlockLduPrecon<vector2> blockVector2Precon;
defineNamedTemplateTypeNameAndDebug(blockVector2Precon, 0);
defineTemplateRunTimeSelectionTable(blockVector2Precon, dictionary);

typedef BlockNoPrecon<vector2> blockNoPreconVector2;
makeBlockPrecon(blockVector2Precon, blockNoPreconVector2);

typedef BlockDiagonalPrecon<vector2> blockDiagonalPreconVector2;
makeBlockPrecon(blockVector2Precon, blockDiagonalPreconVector2);

typedef BlockGaussSeidelPrecon<vector2> blockGaussSeidelPreconVector2;
makeBlockPrecon(blockVector2Precon, blockGaussSeidelPreconVector2);

typedef BlockCholeskyPrecon<vector2> blockCholeskyPreconVector2;
makeBlockPrecon(blockVector2Precon, blockCholeskyPreconVector2);


// Smoothers
typedef BlockLduSmoother<vector2> blockVector2Smoother;
defineNamedTemplateTypeNameAndDebug(blockVector2Smoother, 0);
defineTemplateRunTimeSelectionTable(blockVector2Smoother, dictionary);

typedef BlockGaussSeidelSmoother<vector2> blockGaussSeidelSmootherVector2;
makeBlockSmoother(blockVector2Smoother, blockGaussSeidelSmootherVector2);

typedef BlockILUSmoother<vector2> blockILUSmootherVector2;
makeBlockSmoother(blockVector2Smoother, blockILUSmootherVector2);


// Solvers
typedef BlockLduSolver<vector2> blockVector2Solver;
defineNamedTemplateTypeNameAndDebug(blockVector2Solver, 0);
defineTemplateRunTimeSelectionTable
(
    blockVector2Solver,
    symMatrix
);

defineTemplateRunTimeSelectionTable
(
    blockVector2Solver,
    asymMatrix
);

typedef BlockDiagonalSolver<vector2> blockDiagonalSolverVector2;
defineNamedTemplateTypeNameAndDebug(blockDiagonalSolverVector2, 0);

typedef BlockBiCGStabSolver<vector2> blockBiCGStabSolverVector2;
makeBlockSolverTypeName(blockBiCGStabSolverVector2);
addSolverToBlockMatrix(Vector2, blockBiCGStabSolverVector2, symMatrix);
addSolverToBlockMatrix(Vector2, blockBiCGStabSolverVector2, asymMatrix);

typedef BlockCGSolver<vector2> blockCGSolverVector2;
makeBlockSolverTypeName(blockCGSolverVector2);
addSolverToBlockMatrix(Vector2, blockCGSolverVector2, symMatrix);

typedef BlockGaussSeidelSolver<vector2> blockGaussSeidelSolverVector2;
makeBlockSolverTypeName(blockGaussSeidelSolverVector2);
addSolverToBlockMatrix(Vector2, blockGaussSeidelSolverVector2, symMatrix);

typedef BlockGMRESSolver<vector2> blockGMRESSolverVector2;
makeBlockSolverTypeName(blockGMRESSolverVector2);
addSolverToBlockMatrix(Vector2, blockGMRESSolverVector2, symMatrix);
addSolverToBlockMatrix(Vector2, blockGMRESSolverVector2, asymMatrix);



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
