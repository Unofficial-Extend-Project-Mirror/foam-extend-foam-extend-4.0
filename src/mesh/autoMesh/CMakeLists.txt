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

set(autoHexMesh autoHexMesh)
set(autoHexMeshDriver ${autoHexMesh}/autoHexMeshDriver)

list(APPEND SOURCES
  ${autoHexMeshDriver}/autoLayerDriver.C
  ${autoHexMeshDriver}/autoLayerDriverShrink.C
  ${autoHexMeshDriver}/autoSnapDriver.C
  ${autoHexMeshDriver}/autoRefineDriver.C
  ${autoHexMeshDriver}/autoHexMeshDriver.C
  ${autoHexMeshDriver}/layerParameters/layerParameters.C
  ${autoHexMeshDriver}/refinementParameters/refinementParameters.C
  ${autoHexMeshDriver}/snapParameters/snapParameters.C
  ${autoHexMeshDriver}/pointData/pointData.C
  ${autoHexMesh}/meshRefinement/meshRefinementBaffles.C
  ${autoHexMesh}/meshRefinement/meshRefinement.C
  ${autoHexMesh}/meshRefinement/meshRefinementMerge.C
  ${autoHexMesh}/meshRefinement/meshRefinementProblemCells.C
  ${autoHexMesh}/meshRefinement/meshRefinementRefine.C
  ${autoHexMesh}/refinementSurfaces/refinementSurfaces.C
  ${autoHexMesh}/shellSurfaces/shellSurfaces.C
  ${autoHexMesh}/trackedParticle/trackedParticle.C
  ${autoHexMesh}/trackedParticle/trackedParticleCloud.C
)

add_foam_library(autoMesh SHARED ${SOURCES})

target_link_libraries(autoMesh PUBLIC dynamicFvMesh edgeMesh)
