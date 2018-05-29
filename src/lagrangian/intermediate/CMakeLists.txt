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

set(PARCELS parcels)
set(BASEPARCELS ${PARCELS}/baseClasses)
set(DERIVEDPARCELS ${PARCELS}/derived)

set(CLOUDS clouds)
set(BASECLOUDS ${CLOUDS}/baseClasses)
set(DERIVEDCLOUDS ${CLOUDS}/derived)


#  Parcels
list(APPEND SOURCES
  ${BASEPARCELS}/reactingParcel/reactingParcel.C
)


#  Cloud base classes
list(APPEND SOURCES
  ${BASECLOUDS}/kinematicCloud/kinematicCloud.C
  ${BASECLOUDS}/thermoCloud/thermoCloud.C
  ${BASECLOUDS}/reactingCloud/reactingCloud.C
  ${BASECLOUDS}/reactingMultiphaseCloud/reactingMultiphaseCloud.C
)


#  kinematic parcel sub-models
set(KINEMATICPARCEL ${DERIVEDPARCELS}/basicKinematicParcel)
list(APPEND SOURCES
  ${KINEMATICPARCEL}/basicKinematicParcel.C
  ${KINEMATICPARCEL}/defineBasicKinematicParcel.C
  ${KINEMATICPARCEL}/makeBasicKinematicParcelSubmodels.C
)


#  thermo parcel sub-models
set(THERMOPARCEL ${DERIVEDPARCELS}/basicThermoParcel)
list(APPEND SOURCES
  ${THERMOPARCEL}/basicThermoParcel.C
  ${THERMOPARCEL}/defineBasicThermoParcel.C
  ${THERMOPARCEL}/makeBasicThermoParcelSubmodels.C
)


#  reacting parcel sub-models
set(REACTINGPARCEL ${DERIVEDPARCELS}/BasicReactingParcel)
list(APPEND SOURCES
  ${REACTINGPARCEL}/defineBasicReactingParcel.C
  ${REACTINGPARCEL}/makeBasicReactingParcelSubmodels.C
)


#  reacting multiphase parcel sub-models
set(REACTINGMPPARCEL ${DERIVEDPARCELS}/BasicReactingMultiphaseParcel)
list(APPEND SOURCES
  ${REACTINGMPPARCEL}/defineBasicReactingMultiphaseParcel.C
  ${REACTINGMPPARCEL}/makeBasicReactingMultiphaseParcelSubmodels.C
)


#  bolt-on models
set(RADIATION submodels/addOns/radiation)
list(APPEND SOURCES
  ${RADIATION}/absorptionEmission/cloudAbsorptionEmission/cloudAbsorptionEmission.C
  ${RADIATION}/scatter/cloudScatter/cloudScatter.C
  submodels/Kinematic/PatchInteractionModel/LocalInteraction/patchInteractionData.C
)

set(KINEMATICINJECTION submodels/Kinematic/InjectionModel)
list(APPEND SOURCES
  ${KINEMATICINJECTION}/KinematicLookupTableInjection/kinematicParcelInjectionData.C
  ${KINEMATICINJECTION}/KinematicLookupTableInjection/kinematicParcelInjectionDataIO.C
  ${KINEMATICINJECTION}/KinematicLookupTableInjection/kinematicParcelInjectionDataIOList.C
)

set(THERMOINJECTION submodels/Thermodynamic/InjectionModel)
list(APPEND SOURCES
  ${THERMOINJECTION}/ThermoLookupTableInjection/thermoParcelInjectionData.C
  ${THERMOINJECTION}/ThermoLookupTableInjection/thermoParcelInjectionDataIO.C
  ${THERMOINJECTION}/ThermoLookupTableInjection/thermoParcelInjectionDataIOList.C
)

set(REACTINGINJECTION submodels/Reacting/InjectionModel)
list(APPEND SOURCES
  ${REACTINGINJECTION}/ReactingLookupTableInjection/reactingParcelInjectionData.C
  ${REACTINGINJECTION}/ReactingLookupTableInjection/reactingParcelInjectionDataIO.C
  ${REACTINGINJECTION}/ReactingLookupTableInjection/reactingParcelInjectionDataIOList.C
)

set(REACTINGMPINJECTION submodels/ReactingMultiphase/InjectionModel)
list(APPEND SOURCES
  ${REACTINGMPINJECTION}/ReactingMultiphaseLookupTableInjection/reactingMultiphaseParcelInjectionData.C
  ${REACTINGMPINJECTION}/ReactingMultiphaseLookupTableInjection/reactingMultiphaseParcelInjectionDataIO.C
  ${REACTINGMPINJECTION}/ReactingMultiphaseLookupTableInjection/reactingMultiphaseParcelInjectionDataIOList.C
)


#  data entries
list(APPEND SOURCES
  submodels/IO/DataEntry/makeDataEntries.C
  submodels/IO/DataEntry/polynomial/polynomial.C
  submodels/IO/DataEntry/polynomial/polynomialIO.C
)


#  integration schemes
list(APPEND SOURCES
  IntegrationScheme/makeIntegrationSchemes.C
)


#  particle forces
list(APPEND SOURCES
  particleForces/particleForces.C
)


#  phase properties
list(APPEND SOURCES
  phaseProperties/phaseProperties/phaseProperties.C
  phaseProperties/phaseProperties/phasePropertiesIO.C
  phaseProperties/phasePropertiesList/phasePropertiesList.C
)

add_foam_library(lagrangianIntermediate SHARED ${SOURCES})

target_link_libraries(lagrangianIntermediate PUBLIC reactionThermophysicalModels radiation compressibleRASModels compressibleLESModels pdf liquidMixture solidMixture)
