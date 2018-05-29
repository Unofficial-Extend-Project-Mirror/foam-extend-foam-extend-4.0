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

add_subdirectory(interDyMFoam)
add_subdirectory(interPhaseChangeFoam)
add_subdirectory(multiphaseInterFoam)
add_subdirectory(twoLiquidMixingFoam)
add_subdirectory(interMixingFoam)
add_subdirectory(bubbleFoam)
add_subdirectory(compressibleInterFoam)
add_subdirectory(twoPhaseEulerFoam)
add_subdirectory(MRFInterFoam)
add_subdirectory(barotropicCavitatingFoam)
add_subdirectory(compressibleInterDyMFoam)
add_subdirectory(interFoam)
add_subdirectory(settlingFoam)
add_subdirectory(porousInterFoam)
add_subdirectory(cavitatingFoam)
