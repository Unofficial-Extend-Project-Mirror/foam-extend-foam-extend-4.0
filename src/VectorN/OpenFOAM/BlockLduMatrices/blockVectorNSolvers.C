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
                                                                              \
typedef BlockGMRESSolver<type > block##Type##GMRESSolver;                     \
makeBlockSolverTypeName(block##Type##GMRESSolver);                            \
addSolverToBlockMatrix(Type, block##Type##GMRESSolver, symMatrix);            \
addSolverToBlockMatrix(Type, block##Type##GMRESSolver, asymMatrix);


forAllVectorNTypes(makeSolver)

#undef makeSolver

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
