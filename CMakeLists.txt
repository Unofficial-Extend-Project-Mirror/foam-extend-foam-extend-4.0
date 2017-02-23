# /*-------------------------------------------------------------------------*\#   ========                 |
#   \\      /  F ield         | foam-extend: Open Source CFD
#    \\    /   O peration     | Version:     3.2
#     \\  /    A nd           | Web:         http://www.foam-extend.org
#      \\/     M anipulation  | For copyright notice see file Copyright
# -----------------------------------------------------------------------------
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
#     CMakeLists.txt file for implementing a test harness for the compilation
#     and test of foam-extend-3.2 using Kitware CTest/CMake/CDash
#
#     The results will be submitted to the CDash server identified by the file
#     CTestConfig.cmake
#
# Author
#     Martin Beaudoin, Hydro-Quebec, 2010. All rights reserved
#
#
# \*-------------------------------------------------------------------------*/

cmake_minimum_required (VERSION 2.8)

project(foam-extend C CXX)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

#
# Determine project version & Set-up automatic update during build
#

include(getGitVersion)

add_custom_target(getGitVersion ALL
  COMMAND ${CMAKE_COMMAND} -P cmake/getGitVersion.cmake
  COMMENT "Getting GIT version"
)

set(PROJECT_VERSION ${FOAM_VERSION})


#
# Include helper functions
#

include(FOAMMacros)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)


#
# External dependencies
#

find_package(MPI REQUIRED)
add_library(mpi SHARED IMPORTED)
set_property(TARGET mpi PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${MPI_C_INCLUDE_PATH})
set_property(TARGET mpi PROPERTY IMPORTED_LOCATION ${MPI_LIBRARY})

find_package(ZLIB REQUIRED)
  
find_package(FLEX REQUIRED)

find_package(Git REQUIRED)


#
# Recurse into the source
#

# Set variabable to indicate that we are compiling internally. Used in applications.
set(FOAM_FOUND 1)

# Write something into FOAMTargets.cmake so that we can append to the file.
file(WRITE ${CMAKE_BINARY_DIR}/cmake/FOAMTargets.cmake "#" )

add_subdirectory(src)
add_subdirectory(applications)


#
# Flags and definitions
#

# Define and process the definitions
set(FOAM_PRECISION "double" CACHE STRING "Numerical precision")
set_property(CACHE FOAM_PRECISION PROPERTY STRINGS single double) 
if(FOAM_PRECISION EQUAL "single")
  target_compile_definitions(OSspecific PUBLIC -DWM_SP)
else()
  target_compile_definitions(OSspecific PUBLIC -DWM_SP)
endif()

target_compile_definitions(OSspecific PUBLIC -DNoRepository)

#option(FOAM_PRECISION "help string describing option" -WM_DP)


#
# Exports and install
#

include(GenerateExportHeader)
set(ConfigPackageLocation lib/cmake/foam)

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/cmake/FOAMConfigVersion.cmake"
  VERSION ${FOAM_VERSION}
  COMPATIBILITY AnyNewerVersion
)

configure_package_config_file(cmake/FOAMConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/cmake/FOAMConfig.cmake
  INSTALL_DESTINATION ${ConfigPackageLocation}
)

install(EXPORT FOAMTargets
  FILE FOAMTargets.cmake
  DESTINATION ${ConfigPackageLocation}
)

install(FILES
    cmake/FOAMConfig.cmake
    "${CMAKE_CURRENT_BINARY_DIR}/cmake/FOAMConfigVersion.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake/FOAMMacros.cmake"
  DESTINATION ${ConfigPackageLocation}
  COMPONENT Devel
)

file(GLOB_RECURSE files "${CMAKE_CURRENT_SOURCE_DIR}/etc/*")
install(FILES ${files} DESTINATION etc)

#
# Register with CMake user package registry
#

export(PACKAGE FOAM)


#
# Make a package
#

set(CPACK_DEBIAN_PACKAGE_MAINTAINER "Henrik Rusche")
set(CPACK_PACKAGE_CONTACT "h.rusche@wikki-gmbh.de")
set(CPACK_GENERATOR "DEB")
set(CPACK_PACKAGING_INSTALL_PREFIX "/opt/foam-extend-4.0")
set(CPACK_SOURCE_STRIP_FILES "1")
include(CPack)

