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
  mixtures/basicMixture/basicMixture.C
  mixtures/basicMixture/basicMixtures.C
  basicThermo/basicThermo.C
  psiThermo/basicPsiThermo/basicPsiThermo.C
  psiThermo/basicPsiThermo/newBasicPsiThermo.C
  psiThermo/hPsiThermo/hPsiThermos.C
  psiThermo/hsPsiThermo/hsPsiThermos.C
  psiThermo/ePsiThermo/ePsiThermos.C
  psiThermo/realGasHThermo/realGasHThermos.C
  psiThermo/realGasEThermo/realGasEThermos.C
  rhoThermo/basicRhoThermo/basicRhoThermo.C
  rhoThermo/basicRhoThermo/newBasicRhoThermo.C
  rhoThermo/hRhoThermo/hRhoThermos.C
  rhoThermo/hsRhoThermo/hsRhoThermos.C
  IAPWS_Waterproperties/IAPWSThermo/IAPWS-IF97.C
  IAPWS_Waterproperties/IAPWSThermo/IAPWSThermos.C
  derivedFvPatchFields/fixedEnthalpy/fixedEnthalpyFvPatchScalarField.C
  derivedFvPatchFields/gradientEnthalpy/gradientEnthalpyFvPatchScalarField.C
  derivedFvPatchFields/mixedEnthalpy/mixedEnthalpyFvPatchScalarField.C
  derivedFvPatchFields/fixedInternalEnergy/fixedInternalEnergyFvPatchScalarField.C
  derivedFvPatchFields/gradientInternalEnergy/gradientInternalEnergyFvPatchScalarField.C
  derivedFvPatchFields/mixedInternalEnergy/mixedInternalEnergyFvPatchScalarField.C
  derivedFvPatchFields/ggiEnthalpyJump/ggiEnthalpyJumpFvPatchFields.C
  derivedFvPatchFields/mixingPlaneEnthalpyJump/mixingPlaneEnthalpyJumpFvPatchFields.C
  derivedFvPatchFields/overlapGgiEnthalpyJump/overlapGgiEnthalpyJumpFvPatchFields.C
  derivedFvPatchFields/wallHeatTransfer/wallHeatTransferFvPatchScalarField.C
  derivedFvPatchFields/temperatureDirectedInletOutletVelocity/temperatureDirectedInletOutletVelocityFvPatchVectorField.C
  derivedFvPatchFields/isentropicTotalPressure/isentropicTotalPressureFvPatchScalarField.C
  derivedFvPatchFields/isentropicTotalTemperature/isentropicTotalTemperatureFvPatchScalarField.C
)

add_foam_library(basicThermophysicalModels SHARED ${SOURCES})

target_link_libraries(basicThermophysicalModels PUBLIC finiteVolume specie thermophysicalFunctions)
