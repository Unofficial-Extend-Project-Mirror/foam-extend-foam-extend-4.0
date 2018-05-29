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
  chemistryReaders/chemkinReader/chemkinReader.C
)
FLEX_TARGET(chemkinLexer chemistryReaders/chemkinReader/chemkinLexer.L
  ${CMAKE_CURRENT_BINARY_DIR}/chemkinLexer.C COMPILE_FLAGS "-+")
list(APPEND SOURCES ${FLEX_chemkinLexer_OUTPUTS})
list(APPEND SOURCES
  chemistryReaders/chemistryReader/makeChemistryReaders.C
  mixtures/basicMultiComponentMixture/basicMultiComponentMixture.C
  combustionThermo/hCombustionThermo/hCombustionThermo.C
  combustionThermo/hCombustionThermo/newhCombustionThermo.C
  combustionThermo/hCombustionThermo/hCombustionThermos.C
  combustionThermo/hsCombustionThermo/hsCombustionThermo.C
  combustionThermo/hsCombustionThermo/newhsCombustionThermo.C
  combustionThermo/hsCombustionThermo/hsCombustionThermos.C
  combustionThermo/hhuCombustionThermo/hhuCombustionThermo.C
  combustionThermo/hhuCombustionThermo/newhhuCombustionThermo.C
  combustionThermo/hhuCombustionThermo/hhuCombustionThermos.C
  reactionThermo/hReactionThermo/hReactionThermo.C
  reactionThermo/hReactionThermo/newhReactionThermo.C
  reactionThermo/hReactionThermo/hReactionThermos.C
  reactionThermo/hsReactionThermo/hsReactionThermo.C
  reactionThermo/hsReactionThermo/newhsReactionThermo.C
  reactionThermo/hsReactionThermo/hsReactionThermos.C
  derivedFvPatchFields/fixedUnburntEnthalpy/fixedUnburntEnthalpyFvPatchScalarField.C
  derivedFvPatchFields/gradientUnburntEnthalpy/gradientUnburntEnthalpyFvPatchScalarField.C
  derivedFvPatchFields/mixedUnburntEnthalpy/mixedUnburntEnthalpyFvPatchScalarField.C
)

add_foam_library(reactionThermophysicalModels SHARED ${SOURCES})

target_link_libraries(reactionThermophysicalModels PUBLIC basicThermophysicalModels)
