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
  cellClassification/cellClassification.C
  cellClassification/cellInfo.C
  cellQuality/cellQuality.C
  cellDist/cellDistFuncs.C
  cellDist/patchWave/patchWave.C
  cellDist/wallPoint/wallPoint.C
  cellFeatures/cellFeatures.C
  edgeFaceCirculator/edgeFaceCirculator.C
  polyMeshZipUpCells/polyMeshZipUpCells.C
  primitiveMeshGeometry/primitiveMeshGeometry.C
  meshSearch/meshSearch.C
  PointEdgeWave/PointEdgeWaveName.C
  PointEdgeWave/pointEdgePoint.C
  regionSplit/regionSplit.C
)

set(searchableSurface searchableSurface)
list(APPEND SOURCES
  ${searchableSurface}/distributedTriSurfaceMesh.C
  ${searchableSurface}/searchableBox.C
  ${searchableSurface}/searchableCylinder.C
  ${searchableSurface}/searchablePlane.C
  ${searchableSurface}/searchablePlate.C
  ${searchableSurface}/searchableSphere.C
  ${searchableSurface}/searchableSurface.C
  ${searchableSurface}/searchableSurfaceCollection.C
  ${searchableSurface}/searchableSurfaces.C
  ${searchableSurface}/searchableSurfacesQueries.C
  ${searchableSurface}/searchableSurfaceWithGaps.C
  ${searchableSurface}/triSurfaceMesh.C
)

set(topoSets sets/topoSets)
list(APPEND SOURCES
  ${topoSets}/cellSet.C
  ${topoSets}/topoSet.C
  ${topoSets}/faceSet.C
  ${topoSets}/pointSet.C
  sets/topoSetSource/topoSetSource.C
)

set(cellSources sets/cellSources)
list(APPEND SOURCES
  ${cellSources}/faceToCell/faceToCell.C
  ${cellSources}/fieldToCell/fieldToCell.C
  ${cellSources}/pointToCell/pointToCell.C
  ${cellSources}/shapeToCell/shapeToCell.C
  ${cellSources}/boxToCell/boxToCell.C
  ${cellSources}/regionToCell/regionToCell.C
  ${cellSources}/rotatedBoxToCell/rotatedBoxToCell.C
  ${cellSources}/labelToCell/labelToCell.C
  ${cellSources}/surfaceToCell/surfaceToCell.C
  ${cellSources}/cellToCell/cellToCell.C
  ${cellSources}/nearestToCell/nearestToCell.C
  ${cellSources}/nbrToCell/nbrToCell.C
  ${cellSources}/zoneToCell/zoneToCell.C
  ${cellSources}/setToCell/setToCell.C
  ${cellSources}/sphereToCell/sphereToCell.C
  ${cellSources}/cylinderToCell/cylinderToCell.C
  ${cellSources}/faceZoneToCell/faceZoneToCell.C
)

set(faceSources sets/faceSources)
list(APPEND SOURCES
  ${faceSources}/faceToFace/faceToFace.C
  ${faceSources}/labelToFace/labelToFace.C
  ${faceSources}/cellToFace/cellToFace.C
  ${faceSources}/normalToFace/normalToFace.C
  ${faceSources}/pointToFace/pointToFace.C
  ${faceSources}/patchToFace/patchToFace.C
  ${faceSources}/boundaryToFace/boundaryToFace.C
  ${faceSources}/zoneToFace/zoneToFace.C
  ${faceSources}/setToFace/setToFace.C
  ${faceSources}/boxToFace/boxToFace.C
  ${faceSources}/rotatedBoxToFace/rotatedBoxToFace.C
)

set(pointSources sets/pointSources)
list(APPEND SOURCES
  ${pointSources}/labelToPoint/labelToPoint.C
  ${pointSources}/pointToPoint/pointToPoint.C
  ${pointSources}/cellToPoint/cellToPoint.C
  ${pointSources}/faceToPoint/faceToPoint.C
  ${pointSources}/boxToPoint/boxToPoint.C
  ${pointSources}/surfaceToPoint/surfaceToPoint.C
  ${pointSources}/zoneToPoint/zoneToPoint.C
  ${pointSources}/setToPoint/setToPoint.C
  ${pointSources}/nearestToPoint/nearestToPoint.C
  surfaceSets/surfaceSets.C
)

set(triSurfaceInterfaces triSurface/triSurface/interfaces)
list(APPEND SOURCES
  triSurface/triSurface/triSurface.C
  triSurface/triSurface/triSurfaceAddressing.C
  triSurface/triSurface/stitchTriangles.C
  triSurface/triSurface/treeDataTriSurface/treeDataTriSurface.C
  triSurface/faceTriangulation/faceTriangulation.C
  triSurface/meshTriangulation/meshTriangulation.C
  triSurface/triSurface/geometricSurfacePatch/geometricSurfacePatch.C
  triSurface/triSurface/surfacePatch/surfacePatch.C
  triSurface/triSurface/surfacePatch/surfacePatchIOList.C
  triSurface/tools/labelledTri/sortLabelledTri.C
  triSurface/triSurfaceFields/triSurfaceFields.C
  triSurface/triSurfaceFields/triSurfacePointFields/triSurfacePointFields.C
  triSurface/orientedSurface/orientedSurface.C
  triSurface/surfaceLocation/surfaceLocation.C
  ${triSurfaceInterfaces}/STL/writeSTL.C
  ${triSurfaceInterfaces}/STL/readSTL.C
)
FLEX_TARGET(readSTLASCII ${triSurfaceInterfaces}/STL/readSTLASCII.L
  ${CMAKE_CURRENT_BINARY_DIR}/readSTLASCII.C COMPILE_FLAGS "-+")
list(APPEND SOURCES ${FLEX_readSTLASCII_OUTPUTS})
list(APPEND SOURCES
  ${triSurfaceInterfaces}/STL/readSTLBINARY.C
  ${triSurfaceInterfaces}/GTS/writeGTS.C
  ${triSurfaceInterfaces}/GTS/readGTS.C
  ${triSurfaceInterfaces}/OBJ/readOBJ.C
  ${triSurfaceInterfaces}/OBJ/writeOBJ.C
  ${triSurfaceInterfaces}/SMESH/writeSMESH.C
  ${triSurfaceInterfaces}/OFF/readOFF.C
  ${triSurfaceInterfaces}/OFF/writeOFF.C
  ${triSurfaceInterfaces}/TRI/writeTRI.C
  ${triSurfaceInterfaces}/TRI/readTRI.C
  ${triSurfaceInterfaces}/DX/writeDX.C
  ${triSurfaceInterfaces}/AC3D/readAC.C
  ${triSurfaceInterfaces}/AC3D/writeAC.C
  ${triSurfaceInterfaces}/VTK/writeVTK.C
  ${triSurfaceInterfaces}/NAS/readNAS.C
)


set(booleanOps triSurface/booleanOps)

set(surfaceIntersection ${booleanOps}/surfaceIntersection)
list(APPEND SOURCES
  ${surfaceIntersection}/surfaceIntersection.C
  ${surfaceIntersection}/surfaceIntersectionFuncs.C
  ${surfaceIntersection}/edgeIntersections.C
)

set(booleanSurface ${booleanOps}/booleanSurface)
list(APPEND SOURCES
  ${booleanSurface}/booleanSurface.C
)

set(intersectedSurface ${booleanOps}/intersectedSurface)
list(APPEND SOURCES
  ${intersectedSurface}/intersectedSurface.C
  ${intersectedSurface}/edgeSurface.C
  triSurface/triSurfaceSearch/triSurfaceSearch.C
  triSurface/octreeData/octreeDataTriSurface.C
  triSurface/octreeData/octreeDataTriSurfaceTreeLeaf.C
  triSurface/surfaceFeatures/surfaceFeatures.C
  triSurface/triSurfaceTools/triSurfaceTools.C
  triSurface/triSurfaceTools/geompack/geompack.C
  twoDPointCorrector/twoDPointCorrector.C
  directMapped/directMappedPolyPatch/directMappedPatchBase.C
  directMapped/directMappedPolyPatch/directMappedPolyPatch.C
  directMapped/directMappedPolyPatch/directMappedWallPolyPatch.C
  directMapped/directMappedPointPatch/directMappedPointPatch.C
  directMapped/directMappedPointPatch/directMappedWallPointPatch.C
)

add_foam_library(meshTools SHARED ${SOURCES})

target_link_libraries(meshTools PUBLIC lagrangianBasic decompositionMethods)
