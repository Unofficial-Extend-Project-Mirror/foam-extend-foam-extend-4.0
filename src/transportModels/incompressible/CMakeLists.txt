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
  viscosityModels/viscosityModel/viscosityModel.C
  viscosityModels/viscosityModel/newViscosityModel.C
  viscosityModels/Newtonian/Newtonian.C
  viscosityModels/powerLaw/powerLaw.C
  viscosityModels/CrossPowerLaw/CrossPowerLaw.C
  viscosityModels/BirdCarreau/BirdCarreau.C
  viscosityModels/HerschelBulkley/HerschelBulkley.C
  viscosityModels/freeSurface/freeSurface.C
  transportModel/transportModel.C
  singlePhaseTransportModel/singlePhaseTransportModel.C
  incompressibleTwoPhaseMixture/twoPhaseMixture.C
)

add_foam_library(incompressibleTransportModels SHARED ${SOURCES})

target_link_libraries(incompressibleTransportModels PUBLIC finiteVolume)

target_include_directories(incompressibleTransportModels PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
)
