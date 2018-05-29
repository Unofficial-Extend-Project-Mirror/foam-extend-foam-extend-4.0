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

set(atomicWeights atomicWeights)
set(specie specie)
set(speciesTable speciesTable)
set(equationOfState equationOfState)
set(reactions reaction/reactions)

list(APPEND SOURCES
  ${atomicWeights}/atomicWeights.C
  ${specie}/specie.C
  ${speciesTable}/speciesTable.C
  ${equationOfState}/perfectGas/perfectGas.C
  ${equationOfState}/cubicEquationOfState/redlichKwong/redlichKwong.C
  ${equationOfState}/cubicEquationOfState/aungierRedlichKwong/aungierRedlichKwong.C
  ${equationOfState}/cubicEquationOfState/pengRobinson/pengRobinson.C
  ${equationOfState}/cubicEquationOfState/soaveRedlichKwong/soaveRedlichKwong.C
  ${reactions}/makeChemkinReactions.C
  ${reactions}/makeReactionThermoReactions.C
  ${reactions}/makeLangmuirHinshelwoodReactions.C
)

add_foam_library(specie SHARED ${SOURCES})

target_link_libraries(specie PUBLIC foam)
