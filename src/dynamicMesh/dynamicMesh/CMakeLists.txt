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

set(polyMeshModifiers polyMeshModifiers)

list(APPEND SOURCES
  ${polyMeshModifiers}/polyMeshModifier/polyMeshModifier.C
  ${polyMeshModifiers}/polyMeshModifier/newPolyMeshModifier.C
)

set(attachDetach ${polyMeshModifiers}/attachDetach)
list(APPEND SOURCES
  ${attachDetach}/attachDetach.C
  ${attachDetach}/attachInterface.C
  ${attachDetach}/detachInterface.C
  ${attachDetach}/attachDetachPointMatchMap.C
)

set(layerAdditionRemoval ${polyMeshModifiers}/layerAdditionRemoval)
list(APPEND SOURCES
  ${layerAdditionRemoval}/layerAdditionRemoval.C
  ${layerAdditionRemoval}/setLayerPairing.C
  ${layerAdditionRemoval}/addCellLayer.C
  ${layerAdditionRemoval}/removeCellLayer.C
  ${polyMeshModifiers}/perfectInterface/perfectInterface.C
  ${polyMeshModifiers}/setUpdater/setUpdater.C
)

set(slidingInterface ${polyMeshModifiers}/slidingInterface)
set(enrichedPatch ${slidingInterface}/enrichedPatch)
list(APPEND SOURCES
  ${enrichedPatch}/enrichedPatch.C
  ${enrichedPatch}/enrichedPatchPointMap.C
  ${enrichedPatch}/enrichedPatchFaces.C
  ${enrichedPatch}/enrichedPatchPointPoints.C
  ${enrichedPatch}/enrichedPatchCutFaces.C
  ${enrichedPatch}/enrichedPatchMasterPoints.C
  ${slidingInterface}/slidingInterface.C
  ${slidingInterface}/slidingInterfaceProjectPoints.C
  ${slidingInterface}/coupleSlidingInterface.C
  ${slidingInterface}/slidingInterfaceAttachedAddressing.C
  ${slidingInterface}/slidingInterfaceClearCouple.C
  ${slidingInterface}/decoupleSlidingInterface.C
)

set(repatchCoverage ${polyMeshModifiers}/repatchCoverage)
list(APPEND SOURCES
  ${repatchCoverage}/repatchCoverage.C
)

set(polyhedralRefinement ${polyMeshModifiers}/polyhedralRefinement)
list(APPEND SOURCES
  ${polyhedralRefinement}/polyhedralRefinement.C
  polyTopoChange/polyTopoChange/polyTopoChange.C
  polyTopoChange/polyTopoChange/actions/topoAction/topoActions.C
  polyTopoChange/polyTopoChanger/polyTopoChanger.C
  polyTopoChange/polyTopoChanger/polyTopoChangerChangeMesh.C
  motionSmoother/motionSmoother.C
  motionSmoother/motionSmootherCheck.C
  motionSmoother/polyMeshGeometry/polyMeshGeometry.C
  motionSolver/motionSolver.C
  refinementData/refinementData.C
  refinementData/refinementDistanceData.C
  refinementData/refinementHistory.C
  directTopoChange/directTopoChange/directTopoChange.C
)

set(directActions directTopoChange/directTopoChange/directActions)
list(APPEND SOURCES
  ${directActions}/addPatchCellLayer.C
  ${directActions}/edgeCollapser.C
  ${directActions}/faceCollapser.C
  ${directActions}/hexRef8.C
  ${directActions}/removeCells.C
  ${directActions}/removeFaces.C
  ${directActions}/removePoints.C
  ${directActions}/combineFaces.C
  ${directActions}/localPointRegion.C
  ${directActions}/duplicatePoints.C
)

set(meshCut directTopoChange/meshCut)
list(APPEND SOURCES
  ${meshCut}/cellLooper/cellLooper.C
  ${meshCut}/cellLooper/topoCellLooper.C
  ${meshCut}/cellLooper/geomCellLooper.C
  ${meshCut}/cellLooper/hexCellLooper.C
  ${meshCut}/directions/directions.C
  ${meshCut}/directions/directionInfo/directionInfo.C
  ${meshCut}/edgeVertex/edgeVertex.C
  ${meshCut}/cellCuts/cellCuts.C
  ${meshCut}/splitCell/splitCell.C
  ${meshCut}/refineCell/refineCell.C
  ${meshCut}/wallLayerCells/wallLayerCells.C
  ${meshCut}/wallLayerCells/wallNormalInfo/wallNormalInfo.C
)

set(modifiers ${meshCut}/modifiers)
list(APPEND SOURCES
  ${modifiers}/boundaryCutter/boundaryCutter.C
  ${modifiers}/meshCutter/meshCutter.C
  ${modifiers}/meshCutAndRemove/meshCutAndRemove.C
  ${modifiers}/undoableMeshCutter/undoableMeshCutter.C
  ${modifiers}/refinementIterator/refinementIterator.C
  ${modifiers}/multiDirRefinement/multiDirRefinement.C
  directTopoChange/polyMeshAdder/faceCoupleInfo.C
  directTopoChange/polyMeshAdder/polyMeshAdder.C
  directTopoChange/boundaryMesh/octreeDataFaceList.C
  directTopoChange/boundaryMesh/boundaryPatch.C
  directTopoChange/boundaryMesh/boundaryMesh.C
  directTopoChange/repatchPolyTopoChanger/repatchPolyTopoChanger.C
)

add_foam_library(dynamicMesh SHARED ${SOURCES})

target_link_libraries(dynamicMesh PUBLIC meshTools)
