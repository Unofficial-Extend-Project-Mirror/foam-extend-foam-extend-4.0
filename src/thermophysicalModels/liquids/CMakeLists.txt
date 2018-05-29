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
  liquid/liquid.C
  H2O/H2O.C
  C7H16/C7H16.C
  C12H26/C12H26.C
  C10H22/C10H22.C
  C8H18/C8H18.C
  IC8H18/IC8H18.C
  C4H10O/C4H10O.C
  C2H6O/C2H6O.C
  IDEA/IDEA.C
  aC10H7CH3/aC10H7CH3.C
  bC10H7CH3/bC10H7CH3.C
  C8H10/C8H10.C
  C16H34/C16H34.C
  C9H20/C9H20.C
  C6H6/C6H6.C
  C7H8/C7H8.C
  C6H14/C6H14.C
  C13H28/C13H28.C
  C14H30/C14H30.C
  C3H8/C3H8.C
  C3H6O/C3H6O.C
  C2H6/C2H6.C
  CH3OH/CH3OH.C
  C2H5OH/C2H5OH.C
  Ar/Ar.C
  N2/N2.C
  MB/MB.C
  CH4N2O/CH4N2O.C
  nC3H8O/nC3H8O.C
  iC3H8O/iC3H8O.C
)

add_foam_library(liquids SHARED ${SOURCES})

target_link_libraries(liquids PUBLIC thermophysicalFunctions)
