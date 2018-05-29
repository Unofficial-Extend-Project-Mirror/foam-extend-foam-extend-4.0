# --------------------------------------------------------------------------
#   ========                 |
#   \      /  F ield         | foam-extend: Open Source CFD
#    \    /   O peration     | Version:     4.1
#     \  /    A nd           | Web:         http://www.foam-extend.org
#      \/     M anipulation  | For copyright notice see file Copyright
# --------------------------------------------------------------------------
# License
#     This file is part of foam-extend.
#
#     foam-extend is free software: you can redistribute it and/or modify it
#     under the terms of the GNU General Public License as published by the
#     Free Software Foundation, either version 3 of the License, or (at your
#     option) any later version.
#
#     foam-extend is distributed in the hope that it will be useful, but
#     WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#     General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.
#
# Description
#     CMakeLists.txt file for libraries and applications
#
# Author
#     Henrik Rusche, Wikki GmbH, 2017. All rights reserved
#
#
# --------------------------------------------------------------------------

list(APPEND SOURCES
  coupledLduMatrix/coupledLduMatrix.C
  coupledLduPrecon/coupledLduPrecon/coupledLduPrecon.C
  coupledLduPrecon/noPrecon/coupledNoPrecon.C
  coupledLduPrecon/diagonalPrecon/coupledDiagonalPrecon.C
  coupledLduPrecon/GaussSeidelPrecon/coupledGaussSeidelPrecon.C
  coupledLduPrecon/CholeskyPrecon/coupledCholeskyPrecon.C
  coupledLduSmoother/coupledLduSmoother/coupledLduSmoother.C
  coupledLduSmoother/GaussSeidelSmoother/coupledGaussSeidelSmoother.C
  coupledLduSmoother/iluSmoother/coupledIluSmoother.C
  coupledLduSolver/coupledSolverPerformance/coupledSolverPerformance.C
  coupledLduSolver/coupledLduSolver/coupledLduSolver.C
  coupledLduSolver/iterativeSolver/coupledIterativeSolver.C
  coupledLduSolver/diagonalSolver/coupledDiagonalSolver.C
  coupledLduSolver/smoothSolver/coupledSmoothSolver.C
  coupledLduSolver/cgSolver/coupledCgSolver.C
  coupledLduSolver/bicgSolver/coupledBicgSolver.C
  coupledLduSolver/bicgStabSolver/coupledBicgStabSolver.C
  coupledFvMatrices/coupledFvMatrices.C
  coupledFvMatrices/coupledFvScalarMatrix/coupledFvScalarMatrix.C
)

add_foam_library(coupledLduMatrix SHARED ${SOURCES})

target_link_libraries(coupledLduMatrix PUBLIC finiteVolume)
