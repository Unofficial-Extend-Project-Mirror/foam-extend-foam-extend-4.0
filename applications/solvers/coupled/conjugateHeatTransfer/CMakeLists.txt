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
#     Free Software Foundation, either version 3 of the License, or {at your
#     option} any later version.
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

list(APPEND sources
  thermalModel/thermalModel.C
)
set(thermalLaws thermalModel/thermalLaws)

list(APPEND sources
  ${thermalLaws}/thermalLaw/thermalLaw.C
  ${thermalLaws}/thermalLaw/newThermalLaw.C
  ${thermalLaws}/constantThermal/constantThermal.C
  ${thermalLaws}/multiMaterialThermal/multiMaterialThermal.C
  ${thermalLaws}/multiMaterialZonesThermal/multiMaterialZonesThermal.C
)

set(thermalGaps thermalModel/thermalGaps)

list(APPEND sources
  ${thermalGaps}/thermalGap/thermalGap.C
  ${thermalGaps}/thermalGap/newThermalGap.C
  ${thermalGaps}/constantThermalGap/constantThermalGap.C
  ${thermalGaps}/constantResistanceThermalGap/constantResistanceThermalGap.C
)

set(thermalSource thermalModel/thermalSource)

list(APPEND sources
  ${thermalSource}/thermalSource/thermalSource.C
  ${thermalSource}/thermalSource/newThermalSource.C
  ${thermalSource}/constantThermalSource/constantThermalSource.C
)

set(BCs fvPatchFields)

list(APPEND sources
  ${BCs}/controlledParabolicVelocity/controlledParabolicVelocityFvPatchVectorField.C
  ${BCs}/chtRegionCoupleBase/chtRegionCoupleBase.C
  ${BCs}/chtRcThermalDiffusivity/chtRcThermalDiffusivityFvPatchScalarField.C
  ${BCs}/chtRcTemperature/chtRcTemperatureFvPatchScalarField.C
  ${BCs}/chtRcThermalDiffusivityResistance/chtRcThermalDiffusivityResistanceFvPatchScalarField.C
  ${BCs}/chtRcThermalDiffusivitySlave/chtRcThermalDiffusivitySlaveFvPatchScalarField.C
  ${BCs}/extendedWallHeatTransfer/extendedWallHeatTransferFvPatchScalarField.C
  ${BCs}/externalRadiation/externalRadiationSource/externalRadiationSource.C
  ${BCs}/externalRadiation/externalRadiationSource/newExternalRadiationSource.C
  ${BCs}/externalRadiation/ersConstantFlux/ersConstantFlux.C
  ${BCs}/externalRadiation/ersViewFactor/ersViewFactor.C
  ${BCs}/externalRadiation/ersPlaneToCylinder/ersPlaneToCylinder.C
  ${BCs}/externalRadiation/ersPointSource/ersPointSource.C
  ${BCs}/solidWallHeatFluxTemperature/solidWallHeatFluxTemperatureFvPatchScalarField.C
  ${BCs}/adiabatic/adiabaticFvPatchScalarField.C
)

set(FOs functionObjects)

list(APPEND sources
  ${FOs}/heatFlux/heatFlux/heatFlux.C
  ${FOs}/heatFlux/heatFluxFunctionObject/heatFluxFunctionObject.C
)

add_foam_library(conjugateHeatTransfer SHARED ${sources})

# This is a prototype! Replace second argument with a list of the highest
# dependencies only. Inherited dependencies will be treated automatically.
#target_link_libraries(conjugateHeatTransfer finiteVolume)
