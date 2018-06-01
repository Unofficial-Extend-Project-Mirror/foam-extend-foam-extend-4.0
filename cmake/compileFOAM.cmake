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

#
# Determine project version & Set-up automatic update during build
#

include(getGitVersion)

add_custom_target(getGitVersion ALL
  COMMAND ${CMAKE_COMMAND} -P cmake/getGitVersion.cmake
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
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
set_property(TARGET mpi PROPERTY IMPORTED_LOCATION ${MPI_C_LIBRARIES})
set_property(TARGET mpi PROPERTY INTERFACE_COMPILE_DEFINITIONS OMPI_SKIP_MPICXX)

find_package(ZLIB REQUIRED)

find_package(FLEX REQUIRED)

find_package(Git REQUIRED)

find_package(OpenMP)

# Path to ParaViewConfig.cmake
set(ParaView_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ThirdParty/packages/ParaView-4.4.0/platforms/linux64GccDPOpt/lib/cmake/paraview-4.4)
find_package(ParaView REQUIRED)

find_package(Mesquite REQUIRED)
if(MESQUITE_FOUND)
  add_library(mesquite SHARED IMPORTED)
  set_property(TARGET mesquite PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${MESQUITE_INCLUDE_DIRS})
  set_property(TARGET mesquite PROPERTY IMPORTED_LOCATION ${MESQUITE_LIBRARY})
endif()

find_package(Scotch REQUIRED)
if(SCOTCH_FOUND)
  add_library(scotch SHARED IMPORTED)
  set_property(TARGET scotch PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${SCOTCH_INCLUDE_DIRS})
  set_property(TARGET scotch PROPERTY IMPORTED_LOCATION ${SCOTCH_LIBRARY})
  set_property(TARGET scotch PROPERTY IMPORTED_NO_SONAME TRUE)
endif()

find_package(Metis REQUIRED)
if(METIS_FOUND)
  add_library(metis STATIC IMPORTED)
  set_property(TARGET metis PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${METIS_INCLUDE_DIRS})
  set_property(TARGET metis PROPERTY IMPORTED_LOCATION ${METIS_LIBRARY})
endif()

find_package(ParMetis REQUIRED)
if(PARMETIS_FOUND)
  add_library(parmetis STATIC IMPORTED)
  set_property(TARGET parmetis PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${PARMETIS_INCLUDE_DIRS})
  set_property(TARGET parmetis PROPERTY IMPORTED_LOCATION ${PARMETIS_LIBRARY})
endif()

find_package(ParMGridGen REQUIRED)
if(PARMGRIDGEN_FOUND)
  add_library(parmgridgen SHARED IMPORTED)
  set_property(TARGET parmgridgen PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${PARMGRIDGEN_INCLUDE_DIRS})
  set_property(TARGET parmgridgen PROPERTY IMPORTED_LOCATION ${PARMGRIDGEN_LIBRARY})
  set_property(TARGET parmgridgen PROPERTY IMPORTED_NO_SONAME TRUE)
endif()

find_package(CCMIO REQUIRED)
if(CCMIO_FOUND)
  add_library(ccmio SHARED IMPORTED)
  set_property(TARGET ccmio PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${CCMIO_INCLUDE_DIRS})
  set_property(TARGET ccmio PROPERTY IMPORTED_LOCATION ${CCMIO_LIBRARY})
  set_property(TARGET ccmio PROPERTY INTERFACE_LINK_LIBRARIES ${CCMIO_LIBRARIES})
endif()


#
# Recurse into the source
#

# Set variable to indicate that we are compiling internally. Used in applications.
set(FOAM_FOUND 1)

# Write something into FOAMTargets.cmake so that we can append to the file.
file(WRITE ${CMAKE_BINARY_DIR}/cmake/FOAMTargets.cmake "#" )

add_subdirectory(src)
add_subdirectory(applications)


#
# Set default build type
#

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING
    "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
    FORCE)
endif()


#
# Definitions inherited by all targets
#

# Single/Double precision
set(FOAM_PRECISION "double" CACHE STRING "Numerical precision")
set_property(CACHE FOAM_PRECISION PROPERTY STRINGS single double)
if(FOAM_PRECISION EQUAL "single")
  target_compile_definitions(OSspecific PUBLIC WM_SP)
else()
  target_compile_definitions(OSspecific PUBLIC WM_DP)
endif()

# Label size
set(FOAM_LABEL_SIZE "32" CACHE STRING "Label size")
set_property(CACHE FOAM_LABEL_SIZE PROPERTY STRINGS 32 64)
target_compile_definitions(OSspecific PUBLIC WM_LABEL_SIZE=${FOAM_LABEL_SIZE})

# No Repository
target_compile_definitions(OSspecific PUBLIC NoRepository)

# linux64 - hardcoded for the moment
target_compile_definitions(OSspecific PUBLIC linux64)

# FOAM's full debug mode
if(${CMAKE_BUILD_TYPE} STREQUAL "Debug")
  target_compile_definitions(OSspecific PUBLIC FULLDEBUG)
  target_compile_options(OSspecific PUBLIC -fdefault-inline -ggdb3)
endif()

#target_link_libraries(
#  OSspecific INTERFACE
#  "-Xlinker --add-needed -Xlinker --no-as-needed"
#)

#target_compile_options(
#  OSspecific PUBLIC
#  -Wall -Wextra -Wno-unused-parameter -Wold-style-cast -Wnon-virtual-dtor
#)

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
set(CPACK_PACKAGING_INSTALL_PREFIX "/opt/foam-extend-4.1")
set(CPACK_SOURCE_STRIP_FILES "1")
include(CPack)

