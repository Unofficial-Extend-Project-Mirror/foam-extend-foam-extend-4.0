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

#  RAS turbulence models
list(APPEND SOURCES
  RASModel/RASModel.C
  laminar/laminar.C
  kEpsilon/kEpsilon.C
  RNGkEpsilon/RNGkEpsilon.C
  LaunderSharmaKE/LaunderSharmaKE.C
  LRR/LRR.C
  LaunderGibsonRSTM/LaunderGibsonRSTM.C
  realizableKE/realizableKE.C
  SpalartAllmaras/SpalartAllmaras.C
  kOmegaSST/kOmegaSST.C
)

#  Wall functions
set(wallFunctions derivedFvPatchFields/wallFunctions)

set(alphatWallFunctions ${wallFunctions}/alphatWallFunctions)
list(APPEND SOURCES
  ${alphatWallFunctions}/alphatWallFunction/alphatWallFunctionFvPatchScalarField.C
)

set(mutWallFunctions ${wallFunctions}/mutWallFunctions)
list(APPEND SOURCES
  ${mutWallFunctions}/mutkWallFunction/mutkWallFunctionFvPatchScalarField.C
  ${mutWallFunctions}/mutWallFunction/mutWallFunctionFvPatchScalarField.C
  ${mutWallFunctions}/mutRoughWallFunction/mutRoughWallFunctionFvPatchScalarField.C
  ${mutWallFunctions}/mutSpalartAllmarasWallFunction/mutSpalartAllmarasWallFunctionFvPatchScalarField.C
  ${mutWallFunctions}/mutSpalartAllmarasStandardWallFunction/mutSpalartAllmarasStandardWallFunctionFvPatchScalarField.C
  ${mutWallFunctions}/mutSpalartAllmarasStandardRoughWallFunction/mutSpalartAllmarasStandardRoughWallFunctionFvPatchScalarField.C
  ${mutWallFunctions}/mutLowReWallFunction/mutLowReWallFunctionFvPatchScalarField.C
)

set(epsilonWallFunctions ${wallFunctions}/epsilonWallFunctions)
list(APPEND SOURCES
  ${epsilonWallFunctions}/epsilonWallFunction/epsilonWallFunctionFvPatchScalarField.C
)

set(omegaWallFunctions ${wallFunctions}/omegaWallFunctions)
list(APPEND SOURCES
  ${omegaWallFunctions}/omegaWallFunction/omegaWallFunctionFvPatchScalarField.C
)

set(kqRWallFunctions ${wallFunctions}/kqRWallFunctions)
list(APPEND SOURCES
  ${kqRWallFunctions}/kqRWallFunction/kqRWallFunctionFvPatchFields.C
)

#  Patch fields
list(APPEND SOURCES
  derivedFvPatchFields/turbulentHeatFluxTemperature/turbulentHeatFluxTemperatureFvPatchScalarField.C
  derivedFvPatchFields/turbulentMixingLengthDissipationRateInlet/turbulentMixingLengthDissipationRateInletFvPatchScalarField.C
  derivedFvPatchFields/turbulentMixingLengthFrequencyInlet/turbulentMixingLengthFrequencyInletFvPatchScalarField.C
  derivedFvPatchFields/turbulentTemperatureCoupledBaffle/turbulentTemperatureCoupledBaffleFvPatchScalarField.C
  derivedFvPatchFields/turbulentTemperatureCoupledBaffle/regionProperties.C
  derivedFvPatchFields/turbulentTemperatureCoupledBaffleMixed/turbulentTemperatureCoupledBaffleMixedFvPatchScalarField.C
  backwardsCompatibility/wallFunctions/backwardsCompatibilityWallFunctions.C
)

add_foam_library(compressibleRASModels SHARED ${SOURCES})

target_link_libraries(compressibleRASModels PUBLIC compressibleTurbulenceModel)

target_include_directories(compressibleRASModels PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/derivedFvPatchFields/wallFunctions/mutWallFunctions>
#  $<INSTALL_INTERFACE:include/mylib>
)
