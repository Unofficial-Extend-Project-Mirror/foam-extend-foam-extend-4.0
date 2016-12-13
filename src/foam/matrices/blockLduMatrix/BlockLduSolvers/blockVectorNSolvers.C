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

\*---------------------------------------------------------------------------*/

#include "blockLduMatrices.H"
#include "addToRunTimeSelectionTable.H"

#include "blockLduPrecons.H"
#include "BlockNoPrecon.H"
#include "blockDiagonalPrecons.H"
#include "blockGaussSeidelPrecons.H"
#include "BlockCholeskyPrecon.H"
#include "BlockILUCpPrecon.H"

#include "blockLduSmoothers.H"
#include "blockGaussSeidelSmoothers.H"
#include "BlockILUSmoother.H"
#include "BlockILUCpSmoother.H"

#include "blockLduSolvers.H"
#include "BlockDiagonalSolver.H"
#include "BlockBiCGStabSolver.H"
#include "BlockCGSolver.H"
#include "BlockGaussSeidelSolver.H"
#include "BlockGMRESSolver.H"

// KRJ: 2012-12-15: Multigrid solver
#include "blockAMGSolvers.H"
#include "blockAMGPrecons.H"
#include "blockMatrixCoarsenings.H"
#include "blockMatrixAgglomerations.H"
#include "blockCoeffNorms.H"
#include "blockCoeffTwoNorms.H"
#include "blockCoeffMaxNorms.H"
#include "blockCoeffComponentNorms.H"

#include "VectorTensorNFields.H"
#include "ExpandTensorN.H"
#include "ExpandTensorNField.H"
#include "VectorNFieldTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

#define makeSolver(type, Type, args...)                                       \
/* Preconditioners */                                                         \
typedef BlockLduPrecon<type > block##Type##Precon;                            \
defineNamedTemplateTypeNameAndDebug(block##Type##Precon, 0);                  \
defineTemplateRunTimeSelectionTable(block##Type##Precon, dictionary);         \
                                                                              \
typedef BlockNoPrecon<type > block##Type##NoPrecon;                           \
makeBlockPrecon(block##Type##Precon, block##Type##NoPrecon);                  \
                                                                              \
typedef BlockDiagonalPrecon<type > block##Type##DiagonalPrecon;               \
makeBlockPrecon(block##Type##Precon, block##Type##DiagonalPrecon);            \
                                                                              \
typedef BlockGaussSeidelPrecon<type > block##Type##GaussSeidelPrecon;         \
makeBlockPrecon(block##Type##Precon, block##Type##GaussSeidelPrecon);         \
                                                                              \
typedef BlockCholeskyPrecon<type > block##Type##CholeskyPrecon;               \
makeBlockPrecon(block##Type##Precon, block##Type##CholeskyPrecon);            \
                                                                              \
typedef BlockILUCpPrecon<type > block##Type##ILUCpPrecon;                     \
makeBlockPrecon(block##Type##Precon, block##Type##ILUCpPrecon);               \
                                                                              \
/* Smoothers */                                                               \
typedef BlockLduSmoother<type > block##Type##Smoother;                        \
defineNamedTemplateTypeNameAndDebug(block##Type##Smoother, 0);                \
defineTemplateRunTimeSelectionTable(block##Type##Smoother, dictionary);       \
                                                                              \
typedef BlockGaussSeidelSmoother<type > block##Type##GaussSeidelSmoother;     \
makeBlockSmoother(block##Type##Smoother, block##Type##GaussSeidelSmoother);   \
                                                                              \
typedef BlockILUSmoother<type > block##Type##ILUSmoother;                     \
makeBlockSmoother(block##Type##Smoother, block##Type##ILUSmoother);           \
                                                                              \
typedef BlockILUCpSmoother<type > block##Type##ILUCpSmoother;                 \
makeBlockSmoother(block##Type##Smoother, block##Type##ILUCpSmoother);         \
                                                                              \
                                                                              \
/* Solvers */                                                                 \
typedef BlockLduSolver<type > block##Type##Solver;                            \
defineNamedTemplateTypeNameAndDebug(block##Type##Solver, 0);                  \
defineTemplateRunTimeSelectionTable                                           \
(                                                                             \
    block##Type##Solver,                                                      \
    symMatrix                                                                 \
);                                                                            \
                                                                              \
defineTemplateRunTimeSelectionTable                                           \
(                                                                             \
    block##Type##Solver,                                                      \
    asymMatrix                                                                \
);                                                                            \
                                                                              \
typedef BlockDiagonalSolver<type > block##Type##DiagonalSolver;               \
defineNamedTemplateTypeNameAndDebug(block##Type##DiagonalSolver, 0);          \
                                                                              \
typedef BlockBiCGStabSolver<type > block##Type##BiCGStabSolver;               \
makeBlockSolverTypeName(block##Type##BiCGStabSolver);                         \
addSolverToBlockMatrix(Type, block##Type##BiCGStabSolver, symMatrix);         \
addSolverToBlockMatrix(Type, block##Type##BiCGStabSolver, asymMatrix);        \
                                                                              \
typedef BlockCGSolver<type > block##Type##CGSolver;                           \
makeBlockSolverTypeName(block##Type##CGSolver);                               \
addSolverToBlockMatrix(Type, block##Type##CGSolver, symMatrix);               \
                                                                              \
typedef BlockGaussSeidelSolver<type > block##Type##GaussSeidelSolver;         \
makeBlockSolverTypeName(block##Type##GaussSeidelSolver);                      \
addSolverToBlockMatrix(Type, block##Type##GaussSeidelSolver, symMatrix);      \
addSolverToBlockMatrix(Type, block##Type##GaussSeidelSolver, asymMatrix);     \
                                                                              \
typedef BlockGMRESSolver<type > block##Type##GMRESSolver;                     \
makeBlockSolverTypeName(block##Type##GMRESSolver);                            \
addSolverToBlockMatrix(Type, block##Type##GMRESSolver, symMatrix);            \
addSolverToBlockMatrix(Type, block##Type##GMRESSolver, asymMatrix);           \
                                                                              \
typedef BlockMatrixCoarsening<type > block##Type##MatrixCoarsening;           \
defineNamedTemplateTypeNameAndDebug(block##Type##MatrixCoarsening, 0);        \
defineTemplateRunTimeSelectionTable(block##Type##MatrixCoarsening, matrix);   \
                                                                              \
typedef BlockMatrixAgglomeration<type > block##Type##MatrixAgglomeration;     \
makeBlockMatrixCoarsening(block##Type##MatrixCoarsening, block##Type##MatrixAgglomeration); \
                                                                              \
typedef BlockCoeffNorm<type > block##Type##CoeffNorm;                         \
defineNamedTemplateTypeNameAndDebug(block##Type##CoeffNorm, 0);               \
defineTemplateRunTimeSelectionTable(block##Type##CoeffNorm, dictionary);      \
                                                                              \
typedef BlockCoeffTwoNorm<type > block##Type##CoeffTwoNorm;                   \
makeBlockCoeffNorm(block##Type##CoeffNorm, block##Type##CoeffTwoNorm);        \
                                                                              \
typedef BlockCoeffComponentNorm<type > block##Type##CoeffComponentNorm;       \
makeBlockCoeffNorm(block##Type##CoeffNorm, block##Type##CoeffComponentNorm);  \
                                                                              \
typedef BlockCoeffMaxNorm<type > block##Type##CoeffMaxNorm;                   \
makeBlockCoeffNorm(block##Type##CoeffNorm, block##Type##CoeffMaxNorm);        \
                                                                              \
typedef BlockAMGSolver<type > block##Type##AmgSolver;                         \
makeBlockSolverTypeName(block##Type##AmgSolver);                              \
addSolverToBlockMatrix(Type, block##Type##AmgSolver, symMatrix);              \
addSolverToBlockMatrix(Type, block##Type##AmgSolver, asymMatrix);             \
                                                                              \
typedef BlockAMGPrecon<type > block##Type##AmgPrecon;                         \
makeBlockPrecon(block##Type##Precon, block##Type##AmgPrecon);                 \

forAllVectorNTypes(makeSolver)

#undef makeSolver

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
