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

add_subdirectory(rhoSonicFoam)
add_subdirectory(sonicDyMFoam)
add_subdirectory(steadyCompressibleSRFFoam)
add_subdirectory(realFluidPisoFoam)
add_subdirectory(steadyUniversalMRFFoam)
add_subdirectory(rhoPimpleFoam)
add_subdirectory(sonicFoam)
add_subdirectory(rhoPorousMRFPimpleFoam)
add_subdirectory(rhopSonicFoam)
add_subdirectory(rhoPisoFoam)
add_subdirectory(steadyCompressibleFoam)
add_subdirectory(dbnsFoam)
add_subdirectory(steadyUniversalFoam)
add_subdirectory(steadyCompressibleMRFFoam)
add_subdirectory(rhoSimpleFoam)
add_subdirectory(sonicLiquidFoam)
add_subdirectory(rhoCentralFoam)
add_subdirectory(rhoPorousSimpleFoam)
add_subdirectory(dbnsTurbFoam)
