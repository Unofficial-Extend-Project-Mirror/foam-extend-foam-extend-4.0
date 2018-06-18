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

add_subdirectory(checkSurfaceCurvature)
add_subdirectory(regionCellSets)
add_subdirectory(checkMesh)
add_subdirectory(createBaffles)
add_subdirectory(faceSet)
add_subdirectory(splitMesh)
add_subdirectory(mergeOrSplitBaffles)
add_subdirectory(objToVTK)
add_subdirectory(moveDyMEngineMesh)
add_subdirectory(transformPoints)
add_subdirectory(cellSet)
add_subdirectory(moveEngineMesh)
add_subdirectory(insideCells)
add_subdirectory(subsetMesh)
add_subdirectory(setsToZones)
add_subdirectory(deformedGeom)
add_subdirectory(autoPatch)
add_subdirectory(mergeMeshes)
add_subdirectory(refineMesh)
add_subdirectory(rotateMesh)
add_subdirectory(attachMesh)
add_subdirectory(renumberMesh)
add_subdirectory(splitMeshRegions)
add_subdirectory(deformedTetFemGeom)
add_subdirectory(moveMesh)
add_subdirectory(flattenMesh)
add_subdirectory(zipUpMesh)
add_subdirectory(stitchMesh)
add_subdirectory(tetDecomposition)
add_subdirectory(setSet)
add_subdirectory(mirrorMesh)
add_subdirectory(moveDynamicMesh)
add_subdirectory(pointSet)
add_subdirectory(createPatch)