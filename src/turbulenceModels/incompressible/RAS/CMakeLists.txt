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
  realizableKE/realizableKE.C
  kOmega/kOmega.C
  kOmegaSST/kOmegaSST.C
  SpalartAllmaras/SpalartAllmaras.C
  LRR/LRR.C
  LaunderGibsonRSTM/LaunderGibsonRSTM.C
  LaunderSharmaKE/LaunderSharmaKE.C
  qZeta/qZeta.C
  LienCubicKE/LienCubicKE.C
  LienCubicKELowRe/LienCubicKELowRe.C
  NonlinearKEShih/NonlinearKEShih.C
  LienLeschzinerLowRe/LienLeschzinerLowRe.C
  LamBremhorstKE/LamBremhorstKE.C
  coupledKEpsilon/coupledKEpsilon.C
  coupledKOmegaSST/coupledKOmegaSST.C
)

#  Wall functions
set(wallFunctions derivedFvPatchFields/wallFunctions)

set(nutWallFunctions ${wallFunctions}/nutWallFunctions)
list(APPEND SOURCES
  ${nutWallFunctions}/nutkWallFunction/nutkWallFunctionFvPatchScalarField.C
  ${nutWallFunctions}/nutWallFunction/nutWallFunctionFvPatchScalarField.C
  ${nutWallFunctions}/nutRoughWallFunction/nutRoughWallFunctionFvPatchScalarField.C
  ${nutWallFunctions}/nutSpalartAllmarasWallFunction/nutSpalartAllmarasWallFunctionFvPatchScalarField.C
  ${nutWallFunctions}/nutSpalartAllmarasStandardWallFunction/nutSpalartAllmarasStandardWallFunctionFvPatchScalarField.C
  ${nutWallFunctions}/nutSpalartAllmarasStandardRoughWallFunction/nutSpalartAllmarasStandardRoughWallFunctionFvPatchScalarField.C
  ${nutWallFunctions}/nutLowReWallFunction/nutLowReWallFunctionFvPatchScalarField.C
  ${nutWallFunctions}/nutCWTWallFunction/nutCWTWallFunctionFvPatchScalarField.C
  ${nutWallFunctions}/nutMEWTWallFunction/nutMEWTWallFunctionFvPatchScalarField.C
)

set(epsilonWallFunctions ${wallFunctions}/epsilonWallFunctions)
list(APPEND SOURCES
  ${epsilonWallFunctions}/epsilonWallFunction/epsilonWallFunctionFvPatchScalarField.C
  ${epsilonWallFunctions}/epsilonNonEqWallFunction/epsilonNonEqWallFunctionFvPatchScalarField.C
)

set(omegaWallFunctions ${wallFunctions}/omegaWallFunctions)
list(APPEND SOURCES
  ${omegaWallFunctions}/omegaWallFunction/omegaWallFunctionFvPatchScalarField.C
  ${omegaWallFunctions}/omegaCWTWallFunction/omegaCWTWallFunctionFvPatchScalarField.C
  ${omegaWallFunctions}/omegaMEWTWallFunction/omegaMEWTWallFunctionFvPatchScalarField.C
)

set(kqRWallFunctions ${wallFunctions}/kqRWallFunctions)
list(APPEND SOURCES
  ${kqRWallFunctions}/kqRWallFunction/kqRWallFunctionFvPatchFields.C
  ${kqRWallFunctions}/kNonEqWallFunction/kNonEqWallFunctionFvPatchScalarField.C
)

set(RWallFunctions ${wallFunctions}/RWallFunctions)
list(APPEND SOURCES
  ${RWallFunctions}/RWallFunction/RWallFunctionFvPatchSymmTensorField.C
)

set(wallFunctionFunctionObjects ${wallFunctions}/functionObjects)
list(APPEND SOURCES
  ${wallFunctionFunctionObjects}/pressureGradient/pressureGradientFunctionObject.C
  ${wallFunctionFunctionObjects}/velocityConvection/velocityConvectionFunctionObject.C
)

#  Patch fields
list(APPEND SOURCES
  derivedFvPatchFields/turbulentHeatFluxTemperature/turbulentHeatFluxTemperatureFvPatchScalarField.C
  derivedFvPatchFields/turbulentMixingLengthDissipationRateInlet/turbulentMixingLengthDissipationRateInletFvPatchScalarField.C
  derivedFvPatchFields/turbulentMixingLengthFrequencyInlet/turbulentMixingLengthFrequencyInletFvPatchScalarField.C
  derivedFvPatchFields/nuRatioEpsilon/nuRatioEpsilonFvPatchScalarField.C
  derivedFvPatchFields/nuRatioOmega/nuRatioOmegaFvPatchScalarField.C
  backwardsCompatibility/wallFunctions/backwardsCompatibilityWallFunctions.C
)

add_foam_library(incompressibleRASModels SHARED ${SOURCES})

target_link_libraries(incompressibleRASModels PUBLIC incompressibleTurbulenceModel)

target_include_directories(incompressibleRASModels PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/derivedFvPatchFields/wallFunctions/nutWallFunctions>
#  $<INSTALL_INTERFACE:include/mylib>
)
