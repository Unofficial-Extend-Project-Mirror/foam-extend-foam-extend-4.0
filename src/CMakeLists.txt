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

add_subdirectory(foam)
add_subdirectory(immersedBoundary)
add_subdirectory(overset)
add_subdirectory(lagrangian)
add_subdirectory(edgeMesh)
add_subdirectory(coupledMatrix)
add_subdirectory(turbulenceModels)
add_subdirectory(tetFiniteElement)
add_subdirectory(finiteArea)
add_subdirectory(errorEstimation)
add_subdirectory(randomProcesses)
add_subdirectory(decompositionMethods)
add_subdirectory(conversion)
add_subdirectory(surfMesh)
add_subdirectory(transportModels)
add_subdirectory(multiSolver)
add_subdirectory(OSspecific)
add_subdirectory(equationReader)
add_subdirectory(POD)
add_subdirectory(solidModels)
add_subdirectory(ODE)
add_subdirectory(engine)
add_subdirectory(lduSolvers)
#add_subdirectory(cudaSolvers)
add_subdirectory(mesh)
add_subdirectory(postProcessing)
add_subdirectory(finiteVolume)
add_subdirectory(sampling)
add_subdirectory(fvAgglomerationMethods)
add_subdirectory(meshTools)
add_subdirectory(dynamicMesh)
add_subdirectory(thermophysicalModels)
add_subdirectory(dbns)
add_subdirectory(conjugateHeatTransfer)
