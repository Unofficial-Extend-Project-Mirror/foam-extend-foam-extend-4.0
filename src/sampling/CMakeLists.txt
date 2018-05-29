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
  probes/probes.C
  probes/probesGrouping.C
  probes/probesFunctionObject/probesFunctionObject.C
  sampledSet/coordSet/coordSet.C
  sampledSet/sampledSet/sampledSet.C
  sampledSet/cloud/cloudSet.C
  sampledSet/face/faceOnlySet.C
  sampledSet/curve/curveSet.C
  sampledSet/uniform/uniformSet.C
  sampledSet/midPoint/midPointSet.C
  sampledSet/midPointAndFace/midPointAndFaceSet.C
  sampledSet/sampledSet/sampledSet.C
  sampledSet/sampledSets/sampledSets.C
  sampledSet/sampledSets/sampledSetsGrouping.C
  sampledSet/sampledSetsFunctionObject/sampledSetsFunctionObject.C
)

set(setWriters sampledSet/writers)

list(APPEND SOURCES
  ${setWriters}/writers.C
  ${setWriters}/gnuplot/gnuplotSetWriterRunTime.C
  ${setWriters}/jplot/jplotSetWriterRunTime.C
  ${setWriters}/raw/rawSetWriterRunTime.C
  ${setWriters}/xmgrace/xmgraceSetWriterRunTime.C
  cuttingPlane/cuttingPlane.C
  sampledSurface/sampledPatch/sampledPatch.C
  sampledSurface/sampledPlane/sampledPlane.C
  sampledSurface/isoSurface/isoSurface.C
  sampledSurface/isoSurface/sampledIsoSurface.C
  sampledSurface/isoSurface/isoSurfaceCell.C
  sampledSurface/isoSurface/sampledIsoSurfaceCell.C
  sampledSurface/distanceSurface/distanceSurface.C
  sampledSurface/sampledCuttingPlane/sampledCuttingPlane.C
  sampledSurface/sampledSurface/sampledSurface.C
  sampledSurface/sampledSurfaces/sampledSurfaces.C
  sampledSurface/sampledSurfaces/sampledSurfacesGrouping.C
  sampledSurface/sampledSurfacesFunctionObject/sampledSurfacesFunctionObject.C
  sampledSurface/sampledTriSurfaceMesh/sampledTriSurfaceMesh.C
  sampledSurface/thresholdCellFaces/thresholdCellFaces.C
  sampledSurface/thresholdCellFaces/sampledThresholdCellFaces.C
)

set(surfWriters sampledSurface/writers)

list(APPEND SOURCES
  ${surfWriters}/surfaceWriter.C
  ${surfWriters}/dx/dxSurfaceWriter.C
  ${surfWriters}/foamFile/foamFileSurfaceWriter.C
  ${surfWriters}/nastran/nastranSurfaceWriter.C
  ${surfWriters}/proxy/proxySurfaceWriter.C
  ${surfWriters}/raw/rawSurfaceWriter.C
  ${surfWriters}/starcd/starcdSurfaceWriter.C
  ${surfWriters}/vtk/vtkSurfaceWriter.C
  graphField/writePatchGraph.C
  graphField/writeCellGraph.C
  graphField/makeGraph.C
)

set(meshToMesh meshToMeshInterpolation/meshToMesh)
list(APPEND SOURCES
  ${meshToMesh}/meshToMesh.C
  ${meshToMesh}/calculateMeshToMeshAddressing.C
  ${meshToMesh}/calculateMeshToMeshWeights.C
)

add_foam_library(sampling SHARED ${SOURCES})

target_link_libraries(sampling PUBLIC finiteVolume surfMesh)
