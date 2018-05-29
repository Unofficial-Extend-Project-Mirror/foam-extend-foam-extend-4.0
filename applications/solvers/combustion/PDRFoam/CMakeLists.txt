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
  XiModels/XiModel/XiModel.C
  XiModels/XiModel/newXiModel.C
  XiModels/fixed/fixed.C
  XiModels/algebraic/algebraic.C
  XiModels/transport/transport.C
  XiModels/XiEqModels/XiEqModel/XiEqModel.C
  XiModels/XiEqModels/XiEqModel/newXiEqModel.C
  XiModels/XiEqModels/Gulder/Gulder.C
  XiModels/XiEqModels/instabilityXiEq/instabilityXiEq.C
  XiModels/XiEqModels/SCOPEBlendXiEq/SCOPEBlendXiEq.C
  XiModels/XiEqModels/SCOPEXiEq/SCOPEXiEq.C
  XiModels/XiGModels/XiGModel/XiGModel.C
  XiModels/XiGModels/XiGModel/newXiGModel.C
  XiModels/XiGModels/KTS/KTS.C
  XiModels/XiGModels/instabilityG/instabilityG.C
  PDRModels/turbulence/PDRkEpsilon/PDRkEpsilon.C
  PDRModels/dragModels/PDRDragModel/PDRDragModel.C
  PDRModels/dragModels/PDRDragModel/newPDRDragModel.C
  PDRModels/dragModels/basic/basic.C
  PDRModels/XiEqModels/basicXiSubXiEq/basicXiSubXiEq.C
  PDRModels/XiGModels/basicXiSubG/basicXiSubG.C
  laminarFlameSpeed/SCOPE/SCOPELaminarFlameSpeed.C
  PDRFoam.C
)

# Set minimal environment for external compilation
if(NOT FOAM_FOUND)
  cmake_minimum_required(VERSION 2.8)
  find_package(FOAM REQUIRED)
endif()

add_foam_executable(PDRFoam
  DEPENDS engine compressibleRASModels laminarFlameSpeedModels dynamicFvMesh
  SOURCES ${SOURCES}
)

target_include_directories(PDRFoam PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/XiModels/XiModel>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/XiModels/XiEqModels/XiEqModel>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/XiModels/XiGModels/XiGModel>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/PDRModels/dragModels/PDRDragModel>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/laminarFlameSpeed/SCOPE>
)
