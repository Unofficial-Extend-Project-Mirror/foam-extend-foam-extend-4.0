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
  kineticTheoryModel/kineticTheoryModel.C
  viscosityModel/viscosityModel/kineticTheoryViscosityModel.C
  viscosityModel/viscosityModel/newKineticTheoryViscosityModel.C
  viscosityModel/Gidaspow/GidaspowViscosity.C
  viscosityModel/Syamlal/SyamlalViscosity.C
  viscosityModel/HrenyaSinclair/HrenyaSinclairViscosity.C
  viscosityModel/none/noneViscosity.C
  conductivityModel/conductivityModel/conductivityModel.C
  conductivityModel/conductivityModel/newConductivityModel.C
  conductivityModel/Gidaspow/GidaspowConductivity.C
  conductivityModel/Syamlal/SyamlalConductivity.C
  conductivityModel/HrenyaSinclair/HrenyaSinclairConductivity.C
  radialModel/radialModel/radialModel.C
  radialModel/radialModel/newRadialModel.C
  radialModel/CarnahanStarling/CarnahanStarlingRadial.C
  radialModel/Gidaspow/GidaspowRadial.C
  radialModel/LunSavage/LunSavageRadial.C
  radialModel/SinclairJackson/SinclairJacksonRadial.C
  granularPressureModel/granularPressureModel/granularPressureModel.C
  granularPressureModel/granularPressureModel/newGranularPressureModel.C
  granularPressureModel/Lun/LunPressure.C
  granularPressureModel/SyamlalRogersOBrien/SyamlalRogersOBrienPressure.C
  frictionalStressModel/frictionalStressModel/frictionalStressModel.C
  frictionalStressModel/frictionalStressModel/newFrictionalStressModel.C
  frictionalStressModel/JohnsonJackson/JohnsonJacksonFrictionalStress.C
  frictionalStressModel/Schaeffer/SchaefferFrictionalStress.C
)

add_foam_library(kineticTheoryModel SHARED ${SOURCES})

target_link_libraries(kineticTheoryModel PUBLIC EulerianInterfacialModels)

target_include_directories(kineticTheoryModel PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../interfacialModels/dragModels/dragModel>
)
