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

set(thermophysicalFunction thermophysicalFunction)
set(NSRDSfunctions NSRDSfunctions)
set(APIfunctions APIfunctions)

list(APPEND SOURCES
  ${thermophysicalFunction}/thermophysicalFunction.C
  ${NSRDSfunctions}/NSRDSfunc0/NSRDSfunc0.C
  ${NSRDSfunctions}/NSRDSfunc1/NSRDSfunc1.C
  ${NSRDSfunctions}/NSRDSfunc2/NSRDSfunc2.C
  ${NSRDSfunctions}/NSRDSfunc3/NSRDSfunc3.C
  ${NSRDSfunctions}/NSRDSfunc4/NSRDSfunc4.C
  ${NSRDSfunctions}/NSRDSfunc5/NSRDSfunc5.C
  ${NSRDSfunctions}/NSRDSfunc6/NSRDSfunc6.C
  ${NSRDSfunctions}/NSRDSfunc7/NSRDSfunc7.C
  ${NSRDSfunctions}/NSRDSfunc14/NSRDSfunc14.C
  ${APIfunctions}/APIdiffCoefFunc/APIdiffCoefFunc.C
)

set(freesteam freesteam-2.1)
list(APPEND SOURCES
  ${freesteam}/common.c
  ${freesteam}/steam.c
  ${freesteam}/b23.c
  ${freesteam}/backwards.c
  ${freesteam}/viscosity.c
  ${freesteam}/thcond.c
  ${freesteam}/region1.c
  ${freesteam}/region2.c
  ${freesteam}/region3.c
  ${freesteam}/region4.c
  ${freesteam}/steam_pv.c
  ${freesteam}/steam_ph.c
  ${freesteam}/steam_pT.c
  ${freesteam}/zeroin.c
)

add_foam_library(thermophysicalFunctions SHARED ${SOURCES})

target_link_libraries(thermophysicalFunctions PUBLIC foam)
