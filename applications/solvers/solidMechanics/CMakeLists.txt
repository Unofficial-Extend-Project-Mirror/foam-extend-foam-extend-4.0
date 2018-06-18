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

add_subdirectory(elasticThermalSolidFoam)
add_subdirectory(elasticOrthoAcpSolidFoam)
add_subdirectory(elasticIncrSolidFoam)
add_subdirectory(viscoElasticSolidFoam)
add_subdirectory(elasticSolidFoam)
add_subdirectory(elasticAcpSolidFoam)
add_subdirectory(icoFsiElasticNonLinULSolidFoam)
add_subdirectory(elasticOrthoNonLinULSolidFoam)
add_subdirectory(utilities)
add_subdirectory(deprecatedSolvers)
add_subdirectory(elasticPlasticNonLinTLSolidFoam)
add_subdirectory(elasticNonLinTLSolidFoam)
add_subdirectory(elasticOrthoSolidFoam)
add_subdirectory(stressFemFoam)
add_subdirectory(elasticNonLinULSolidFoam)
add_subdirectory(elasticNonLinIncrTLSolidFoam)
add_subdirectory(elasticPlasticSolidFoam)
add_subdirectory(elasticIncrAcpSolidFoam)
add_subdirectory(elasticPlasticNonLinULSolidFoam)