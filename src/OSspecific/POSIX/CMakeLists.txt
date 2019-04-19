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

get_property(inc TARGET foam PROPERTY INCLUDE_DIRECTORIES)
include_directories(${inc})

list(APPEND SOURCES
  signals/sigFpe.C
  signals/sigSegv.C
  signals/sigInt.C
  signals/sigQuit.C
  regExp.C
  timer.C
  fileStat.C
  POSIX.C
  cpuTime/cpuTime.C
  clockTime/clockTime.C
  memInfo/memInfo.C
  multiThreader/multiThreader.C

# Note: fileMonitor assumes inotify by default. Compile with -DFOAM_USE_STAT
# to use stat (=timestamps) instead of inotify
 fileMonitor.C
)

#ifdef SunOS64
#list(APPEND SOURCES
#  dummyPrintStack.C
#)
#else
list(APPEND SOURCES
  printStack.C
)
#endif

add_foam_library(OSspecific SHARED ${SOURCES})

add_dependencies(OSspecific foam_lnInclude)
