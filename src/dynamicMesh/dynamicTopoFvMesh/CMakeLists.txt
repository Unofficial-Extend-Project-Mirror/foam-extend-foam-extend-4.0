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
  eMesh/eMesh.C
  eMesh/eMeshDemandDrivenData.C
  eMesh/eBoundaryMesh/eBoundaryMesh.C
)

set(ePatches eMesh/ePatches)
list(APPEND SOURCES
  ${ePatches}/ePatch/ePatch.C
  ${ePatches}/ePatch/newEPatch.C
  dynamicTopoFvMesh.C
  dynamicTopoFvMeshCheck.C
  dynamicTopoFvMeshReOrder.C
  dynamicTopoFvMeshMapping.C
  edgeSwap.C
  edgeBisect.C
  edgeCollapse.C
  coupledMesh/coupleMap.C
  coupledMesh/dynamicTopoFvMeshCoupled.C
  coupledMesh/subMeshProcessorPolyPatch.C
  coupledMesh/subMeshProcessorFvPatch.C
  convexSetAlgorithm/convexSetAlgorithm.C
  convexSetAlgorithm/faceSetAlgorithm.C
  convexSetAlgorithm/cellSetAlgorithm.C
  fieldMapping/topoMapper.C
  fieldMapping/fluxCorrector.C
  fieldMapping/topoCellMapper.C
  fieldMapping/topoPatchMapper.C
  fieldMapping/topoSurfaceMapper.C
)

set(tetMetrics tetMetrics)
list(APPEND SOURCES
  ${tetMetrics}/tetMetric.C
  ${tetMetrics}/tetMetrics.C
)

set(lengthScaleEstimator lengthScaleEstimator)
list(APPEND SOURCES
  ${lengthScaleEstimator}/lengthScaleEstimator.C
)

add_foam_library(dynamicTopoFvMesh SHARED ${SOURCES})

target_link_libraries(dynamicTopoFvMesh PUBLIC dynamicFvMesh)
