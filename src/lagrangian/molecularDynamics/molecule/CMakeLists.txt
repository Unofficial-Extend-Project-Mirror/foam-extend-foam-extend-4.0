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

set(interactionLists interactionLists)
set(referredMolecule ${interactionLists}/referredMolecule)
set(referredCellList ${interactionLists}/referredCellList)
set(referredCell ${interactionLists}/referredCell)
set(referralLists ${interactionLists}/referralLists)
set(directInteractionList ${interactionLists}/directInteractionList)

list(APPEND SOURCES
  ${referralLists}/sendingReferralList.C
  ${referralLists}/receivingReferralList.C
  ${referredCellList}/referredCellList.C
  ${referredCell}/referredCell.C
  ${referredMolecule}/referredMolecule.C
  ${directInteractionList}/directInteractionList.C
  ${interactionLists}/interactionLists.C
)

set(reducedUnits reducedUnits)
list(APPEND SOURCES
  ${reducedUnits}/reducedUnits.C
  ${reducedUnits}/reducedUnitsIO.C
)

set(molecule molecule)
list(APPEND SOURCES
  ${molecule}/molecule.C
  ${molecule}/moleculeIO.C
)

set(moleculeCloud moleculeCloud)
list(APPEND SOURCES
  ${moleculeCloud}/moleculeCloud.C
)

add_foam_library(molecule SHARED ${SOURCES})

target_link_libraries(molecule PUBLIC potential)
