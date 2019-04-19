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

add_subdirectory(netgenNeutralToFoam)
add_subdirectory(writeMeshObj)
add_subdirectory(polyDualMesh)
add_subdirectory(gmshToFoam)
add_subdirectory(cfx4ToFoam)
add_subdirectory(ideasUnvToFoam)
add_subdirectory(foamMeshToElmer)
add_subdirectory(fluent3DMeshToElmer)
add_subdirectory(sammToFoam)
add_subdirectory(Optional)
add_subdirectory(mshToFoam)
add_subdirectory(foamMeshToAbaqus)
add_subdirectory(starToFoam)
add_subdirectory(kivaToFoam)
add_subdirectory(foamMeshToFluent)
add_subdirectory(gambitToFoam)
add_subdirectory(ansysToFoam)
add_subdirectory(foamToStarMesh)
add_subdirectory(plot3dToFoam)
add_subdirectory(fluent3DMeshToFoam)
add_subdirectory(tetgenToFoam)
add_subdirectory(fluentMeshToFoam)
add_subdirectory(star4ToFoam)
