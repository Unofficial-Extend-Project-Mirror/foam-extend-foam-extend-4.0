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
  rheologyModel/rheologyModel.C
)
set(rheologyLaws rheologyModel/rheologyLaws)

list(APPEND SOURCES
  ${rheologyLaws}/rheologyLaw/rheologyLaw.C
  ${rheologyLaws}/rheologyLaw/newRheologyLaw.C
  ${rheologyLaws}/linearElastic/linearElastic.C
  ${rheologyLaws}/multiMaterial/multiMaterial.C
  ${rheologyLaws}/MaxwellViscoelastic/MaxwellViscoelastic.C
  ${rheologyLaws}/BurgersViscoelastic/BurgersViscoelastic.C
  ${rheologyLaws}/KelvinSLSViscoelastic/KelvinSLSViscoelastic.C
  ${rheologyLaws}/MaxwellSLSViscoelastic/MaxwellSLSViscoelastic.C
  ${rheologyLaws}/MaxwellElasticViscoelastic/MaxwellElasticViscoelastic.C
  ${rheologyLaws}/PronyViscoelastic/PronyViscoelastic.C
  thermalModel/thermalModel.C
)
set(thermalLaws thermalModel/thermalLaws)

list(APPEND SOURCES
  ${thermalLaws}/thermalLaw/thermalLaw.C
  ${thermalLaws}/thermalLaw/newThermalLaw.C
  ${thermalLaws}/constantThermal/constantThermal.C
  ${thermalLaws}/multiMaterialThermal/multiMaterialThermal.C
  cohesiveLaws/cohesiveLaw/cohesiveLaw.C
  cohesiveLaws/Dugdale/DugdaleCohesiveLaw.C
  cohesiveLaws/linear/linearCohesiveLaw.C
  fvPatchFields/tractionDisplacement/tractionDisplacementFvPatchVectorField.C
  fvPatchFields/tractionDisplacementThermo/tractionDisplacementThermoFvPatchVectorField.C
  fvPatchFields/nusselt/nusseltFvPatchScalarField.C
  fvPatchFields/newDirectionMixed/newDirectionMixedFvPatchFields.C
  fvPatchFields/cohesiveLaw/cohesiveLawFvPatchVectorField.C
  fvPatchFields/cohesiveZone/cohesiveZoneFvPatchVectorField.C
)

add_foam_library(materialModels SHARED ${SOURCES})

target_link_libraries(materialModels PUBLIC finiteVolume)
