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

add_custom_target(ImportMessage
  ${CMAKE_COMMAND} -E cmake_echo_color --red --bold
  "Using FOAM ${FOAM_VERSION} in ${CMAKE_CURRENT_LIST_DIR}"
)

function(add_foam_library lib)
  set(options USERSPACE)
  set(oneValueArgs "")
  set(multiValueArgs "")
  cmake_parse_arguments(AFL "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

  # Create target for lnInclude and use it
  include_directories(lnInclude)
  add_custom_command(
    OUTPUT lnInclude/uptodate
    DEPENDS CMakeLists.txt
    COMMAND rm -rf lnInclude
    COMMAND wmakeLnInclude .
    COMMAND touch lnInclude/uptodate
  )
  add_custom_target(${lib}_lnInclude DEPENDS lnInclude/uptodate)

  if(${AFL_USERSPACE})
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY $ENV{FOAM_USER_LIBBIN})
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY $ENV{FOAM_USER_APPBIN})
  endif()

  # Add the library to the targets and set include paths
  add_library(${lib} ${AFL_UNPARSED_ARGUMENTS})
  add_dependencies(${lib} ${lib}_lnInclude)
  target_include_directories(${lib} PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/lnInclude>
    $<INSTALL_INTERFACE:include/${lib}>
  )

  # Install library and headers to install location
  install(TARGETS ${lib} DESTINATION lib EXPORT FOAMTargets)
  file(GLOB _files FOLLOW_SYMLINKS "${CMAKE_CURRENT_SOURCE_DIR}/lnInclude/*")
  set (_resolvedFiles "")
  foreach (_file ${_files})
    get_filename_component(_resolvedFile "${_file}" REALPATH)
    list (APPEND _resolvedFiles "${_resolvedFile}")
  endforeach()
  install(FILES ${_resolvedFiles} DESTINATION include/${lib})

  # Export target to include them from externals builds to build location
  export(TARGETS ${lib} APPEND FILE ${CMAKE_BINARY_DIR}/cmake/FOAMTargets.cmake)

  #generate_export_header(${lib})
  set_property(TARGET ${lib} PROPERTY VERSION ${FOAM_VERSION})
  set_property(TARGET ${lib} PROPERTY SOVERSION 3)
  set_property(TARGET ${lib} PROPERTY INTERFACE_FOAM_MAJOR_VERSION 4)
  set_property(TARGET ${lib} APPEND PROPERTY
    COMPATIBLE_INTERFACE_STRING FOAM_MAJOR_VERSION
  )

  if(NOT DEFINED HAVE_FOAM)
    add_dependencies(${lib} ImportMessage)
  endif()
endfunction()

function(add_foam_executable exe)
  set(options USERSPACE)
  set(oneValueArgs "")
  set(multiValueArgs DEPENDS SOURCES)
  cmake_parse_arguments(AFE "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

  if(${AFE_USERSPACE})
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY $ENV{FOAM_USER_LIBBIN})
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY $ENV{FOAM_USER_APPBIN})
  endif()

  add_executable(${exe} ${AFE_SOURCES})
  target_link_libraries(${exe}
    "-Xlinker --add-needed -Xlinker --no-as-needed"
  )
  target_link_libraries(${exe} ${AFE_DEPENDS})
  install (TARGETS ${exe} DESTINATION bin)

  if(NOT DEFINED HAVE_FOAM)
    add_dependencies(${exe} ImportMessage)
  endif()
endfunction()
