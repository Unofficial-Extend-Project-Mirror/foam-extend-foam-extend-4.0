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

if(NOT GIT_FOUND)
  find_package(Git QUIET)
endif()

if(GIT_FOUND)
  # Try to get version from from git
  execute_process(
    COMMAND ${GIT_EXECUTABLE} describe --tags --dirty=-dirty
    OUTPUT_VARIABLE GIT_VERSION
    ERROR_VARIABLE dummy
    RESULT_VARIABLE res
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )

  if(res EQUAL 0)
    string(REPLACE "g" "" GIT_VERSION "${GIT_VERSION}")
    string(REGEX REPLACE "^v([0-9]+\\.?[0-9]*\\.?([0-9]*)).*" "\\1" FOAM_VERSION "${GIT_VERSION}")
  endif()
endif()

if(NOT GIT_VERSION)
  if(EXISTS cmake/storedVersion.cmake)
    # Fall-back to stored version
    include(cmake/storedVersion.cmake)
  else()
    # We should not be here. Set some defaults!
    set(GIT_VERSION "4.0-unknown")
    set(FOAM_VERSION "4.0")
  endif()
endif()

# Save version information so that it can be exported with the project
# and included as a fall-back - see above
configure_file(cmake/storedVersion.cmake.in
  cmake/storedVersion.cmake
)

# Configure global.C
configure_file(${CMAKE_SOURCE_DIR}/src/foam/global/global.C.in
  src/foam/global/global.C
)

