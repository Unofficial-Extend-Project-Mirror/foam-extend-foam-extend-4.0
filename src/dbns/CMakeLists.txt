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
  dbnsFlux/roeFlux/roeFlux.C
  dbnsFlux/rusanovFlux/rusanovFlux.C
  dbnsFlux/betaFlux/betaFlux.C
  dbnsFlux/hllcFlux/hllcFlux.C
  dbnsFlux/hllcALEFlux/hllcALEFlux.C
  basicNumericFlux/basicNumericFlux.C
  basicNumericFlux/newBasicNumericFlux.C
  numericFlux/numericFluxes.C
  multigrid/mgMeshLevel/mgMeshLevel.C
  multigrid/mgMeshLevel/fineMgMeshLevel.C
  multigrid/mgMeshLevel/coarseMgMeshLevel.C
  multigrid/mgFieldLevel/mgFieldLevel.C
  multigrid/mgFieldLevel/fineMgFieldLevel.C
  multigrid/mgFieldLevel/coarseMgFieldLevel.C
  timeStepping/localTimeStep/localTimeStep.C
  timeStepping/backwardDualDdtScheme/backwardDualDdtSchemes.C
  timeStepping/EulerLocalDdtScheme/EulerLocalDdtSchemes.C
)

if(PARMGRIDGEN_FOUND)
  add_foam_library(dbns SHARED ${SOURCES})

  target_link_libraries(dbns PUBLIC compressibleTurbulenceModel parmgridgen)
else()
  message(WARNING dbns will not be build)
endif()
