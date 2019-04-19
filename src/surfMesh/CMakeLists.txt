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
  surfZone/surfZone/surfZone.C
  surfZone/surfZone/surfZoneIOList.C
  surfZone/surfZoneIdentifier/surfZoneIdentifier.C
  MeshedSurfaceAllocator/MeshedSurfaceIOAllocator.C
  MeshedSurface/MeshedSurfaceCore.C
  MeshedSurface/MeshedSurfaces.C
  UnsortedMeshedSurface/UnsortedMeshedSurfaces.C
  MeshedSurfaceProxy/MeshedSurfaceProxyCore.C
  surfaceRegistry/surfaceRegistry.C
  surfMesh/surfMesh.C
  surfMesh/surfMeshClear.C
  surfMesh/surfMeshIO.C
  surfFields/surfFields/surfFields.C
  surfFields/surfPointFields/surfPointFields.C
)

set(surfaceFormats surfaceFormats)
list(APPEND SOURCES
  ${surfaceFormats}/surfaceFormatsCore.C
  ${surfaceFormats}/ac3d/AC3DsurfaceFormatCore.C
  ${surfaceFormats}/ac3d/AC3DsurfaceFormatRunTime.C
  ${surfaceFormats}/ftr/FTRsurfaceFormatRunTime.C
  ${surfaceFormats}/gts/GTSsurfaceFormatRunTime.C
  ${surfaceFormats}/nas/NASsurfaceFormatCore.C
  ${surfaceFormats}/nas/NASsurfaceFormatRunTime.C
  ${surfaceFormats}/obj/OBJsurfaceFormatRunTime.C
  ${surfaceFormats}/off/OFFsurfaceFormatRunTime.C
  ${surfaceFormats}/ofs/OFSsurfaceFormatCore.C
  ${surfaceFormats}/ofs/OFSsurfaceFormatRunTime.C
  ${surfaceFormats}/smesh/SMESHsurfaceFormatRunTime.C
  ${surfaceFormats}/starcd/STARCDsurfaceFormatCore.C
  ${surfaceFormats}/starcd/STARCDsurfaceFormatRunTime.C
  ${surfaceFormats}/stl/STLsurfaceFormatCore.C
  ${surfaceFormats}/stl/STLsurfaceFormatRunTime.C
)
FLEX_TARGET(STLsurfaceFormatASCII ${surfaceFormats}/stl/STLsurfaceFormatASCII.L
  ${CMAKE_CURRENT_BINARY_DIR}/STLsurfaceFormatASCII.C COMPILE_FLAGS "-+")
list(APPEND SOURCES ${FLEX_STLsurfaceFormatASCII_OUTPUTS})
list(APPEND SOURCES
  ${surfaceFormats}/tri/TRIsurfaceFormatCore.C
  ${surfaceFormats}/tri/TRIsurfaceFormatRunTime.C
  ${surfaceFormats}/vtk/VTKsurfaceFormatCore.C
  ${surfaceFormats}/vtk/VTKsurfaceFormatRunTime.C
  ${surfaceFormats}/wrl/WRLsurfaceFormatCore.C
  ${surfaceFormats}/wrl/WRLsurfaceFormatRunTime.C
  ${surfaceFormats}/x3d/X3DsurfaceFormatCore.C
  ${surfaceFormats}/x3d/X3DsurfaceFormatRunTime.C
)

add_foam_library(surfMesh SHARED ${SOURCES})

target_link_libraries(surfMesh PUBLIC foam)
