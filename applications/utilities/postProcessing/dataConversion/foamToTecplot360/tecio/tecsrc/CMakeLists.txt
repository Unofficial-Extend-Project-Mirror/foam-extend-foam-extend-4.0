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
  alloc.cpp
  arrlist.cpp
  auxdata.cpp
  dataio4.cpp
  dataio.cpp
  dataset0.cpp
  dataset.cpp
  datautil.cpp
  filestream.cpp
  geom2.cpp
  q_msg.cpp
  q_unicode.cpp
  set.cpp
  strlist.cpp
  strutil.cpp
  tassert.cpp
  tecxxx.cpp
  TranslatedString.cpp
)


add_foam_library(tecio SHARED ${SOURCES})


target_compile_definitions(tecio PRIVATE MAKEARCHIVE USEENUM THREED ENGINE)
if(UNIX)
  target_compile_definitions(tecio PRIVATE LINUX)
  target_compile_options(tecio PRIVATE -U_WIN32)
  if(CMAKE_${LANG}_COMPILER_ID STREQUAL "Intel")
    target_compile_definitions(tecio PRIVATE LINUXI64)
  else()
    target_compile_definitions(tecio PRIVATE LINUX64)
    target_compile_options(tecio PRIVATE -U_WIN32)
  endif()
elseif(APPLE)
  target_compile_definitions(tecio PRIVATE DARWIN LINUX64)
  target_compile_options(tecio PRIVATE -U_WIN32)
elseif(MINGW)
  target_compile_definitions(tecio PRIVATE LINUX)
endif()

# This tries to mimic
#
#if defined(linux64)
#    TECIO_FLAGS = -DMAKEARCHIVE -DLINUX -DLINUX64 -DUSEENUM -DTHREED -U_WIN32 -DENGINE
#elif defined(darwinIntel) || defined(darwinIntel64) || defined(darwinPpc) || defined(darwinPpc64)
#    TECIO_FLAGS = -DMAKEARCHIVE -DDARWIN -DUSEENUM -DTHREED -U_WIN32 -DENGINE
#elif defined(linuxIA64)
#    TECIO_FLAGS = -DMAKEARCHIVE -DLINUX -DLINUXI64 -DUSEENUM -DTHREED -U_WIN32 -DENGINE
#elif defined(linux)
#    TECIO_FLAGS = -DMAKEARCHIVE -DLINUX -DLINUX64 -DUSEENUM -DTHREED -U_WIN32 -DENGINE
#elif defined(darwinIntel)
#   warning "DarwinIntel architecture detected"
#    TECIO_FLAGS = -DMAKEARCHIVE -DLINUX -DLINUX64 -DUSEENUM -DTHREED -U_WIN32 -DENGINE
#elif defined(mingw)
#   warning "MSWin mingw architecture detected"
#    TECIO_FLAGS = -DMAKEARCHIVE -DLINUX -DUSEENUM -DTHREED -DENGINE
#else
#error "architecture not supported for compiling tecio."
#endif
