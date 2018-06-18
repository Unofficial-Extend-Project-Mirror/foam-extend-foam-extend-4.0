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

set(lduPrecon lduPrecon)
list(APPEND SOURCES
  ${lduPrecon}/CholeskyPrecon/CholeskyPrecon.C
  ${lduPrecon}/ILU0/ILU0.C
  ${lduPrecon}/ILUC0/ILUC0.C
  ${lduPrecon}/ILUCp/ILUCp.C
  ${lduPrecon}/symGaussSeidelPrecon/symGaussSeidelPrecon.C
  ${lduPrecon}/amgPrecon/amgPrecon.C
)

set(lduSmoother lduSmoother)
list(APPEND SOURCES
  ${lduSmoother}/symGaussSeidelSmoother/symGaussSeidelSmoother.C
  ${lduSmoother}/iluSmoother/iluSmoother.C
  ${lduSmoother}/iluC0Smoother/iluC0Smoother.C
)

set(lduSolver lduSolver)
list(APPEND SOURCES
  ${lduSolver}/cgSolver/cgSolver.C
  ${lduSolver}/bicgSolver/bicgSolver.C
  ${lduSolver}/bicgStabSolver/bicgStabSolver.C
  ${lduSolver}/gmresSolver/gmresSolver.C
  ${lduSolver}/amgSolver/amgSolver.C
  ${lduSolver}/fpeAmgSolver/fpeAmgSolver.C
  ${lduSolver}/mpeAmgSolver/mpeAmgSolver.C
  ${lduSolver}/rreAmgSolver/rreAmgSolver.C
  ${lduSolver}/deflationSolver/deflationSolver.C
)

set(amg amg)
list(APPEND SOURCES
  ${amg}/amgCycle.C
  ${amg}/fineAmgLevel.C
  ${amg}/coarseAmgLevel.C
)

set(amgPolicy ${amg}/amgPolicy)
list(APPEND SOURCES
  ${amgPolicy}/amgPolicy.C
  ${amgPolicy}/pamgPolicy.C
  ${amgPolicy}/clusterAmgPolicy.C
  ${amgPolicy}/selectionAmgPolicy.C
)

add_foam_library(lduSolvers SHARED ${SOURCES})

target_link_libraries(lduSolvers PUBLIC foam)
