#------------------------------------------------------------------------------
# =========                 |
# \\      /  F ield         | foam-extend: Open Source CFD
#  \\    /   O peration     | Version:     3.2
#   \\  /    A nd           | Web:         http://www.foam-extend.org
#    \\/     M anipulation  | For copyright notice see file Copyright
#------------------------------------------------------------------------------
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
# Script
#     doc/buildInstructions/Windows/README.txt
#
# Description
#     Environment setup and build instructions for MinGW-based Windows builds.
#
# Author:
#     Cesare Guardino, Alstom Power Ltd., (2015)
#
#------------------------------------------------------------------------------

1. INTRODUCTION
===============

It is strongly recommended to install all required systems tools and utilities
in a directory with no embedded white spaces. You can skip some steps if you
already have the correct tools installed on your system.


2. INSTRUCTIONS FOR BUILDING ON WINDOWS
=======================================

1) Download and install 7-Zip (see reference list below). This is necessary to
   be able to unzip the rest of the required packages mentioned in later steps.

2) Download and install wget, CMake, Git, MinGW-w64, ParaView, PExports, Java,
   Subversion, Python (see reference list below). Note that some components
   (for example PExports and GMake) may already be available in other packages
   (such as Strawberry Perl), although you need to be careful of the
    version numbers.

3) Download and install MSYS (see reference list below).  If this is your first
   use of MSYS, create a file fstab in c/MinGW/msys/1.0/etc with single-line
   contents:

   c:/mingw /mingw

   to mount your home directory (this assumes you have installed MSYS to the default
   c:/MinGW directory).
   From Windows, you will find your home directory under
   c:/MinGW/msys/1.0/home/<user name>

4) Download the foam-extend source code from

   http://sourceforge.net/projects/openfoam-extend/

   using the command:

   git clone --branch nextRelease git://git.code.sf.net/p/openfoam-extend/foam-extend-3.1 foam-extend-3.2

5) Open a new MSYS shell and chdir to your foam-extend-3.2 source directory.

6) Edit the user-modifiable entries in etc/bashrc.mingw to point to the
   locations where you have installed the required system tools in the first two
   steps (the first two functions only in bashrc.mingw). The order in which they
   are added to the PATH (in the add_to_path function) is very important. The
   PATH is read from left to right, so ensure OpenMPI and MinGW compiler are
   prepended last to the PATH, so they are found first.

7) Source the file edited in the previous step using the command:
       . etc/bashrc.mingw
   (notice the dot "." command to source a file). This action sources other
   files as required, and also performs a check on the versions of the installed
   system tools.  An example output is shown at the end of this file.

8) If "mpirun.exe" failed to run due to a missing DLL (eg. msvcr100.dll), you
   may need to install the Microsoft Visual C++ 2010 Redistributable Package (see
   reference list below).

9) If the version checks were all successful and printed the expected version
   details, proceed with compiling the code. Run "Allwmake.mingw". This will
   download and build all required third-party dependencies, build the entire
   foam-extend code, and also create a stand-alone .zip package.

10) After the build has completed, you can run foam in either of two ways:

   (a). From the MSYS shell. This allows use of the utility programs and shell
        scripts in the bin directory (such as paraFoam). This is the recommended
        approach for developers.

   (b). From a standard Windows CMD.exe command prompt using the created
        foam-extend-3.2-win-x64.zip stand-alone package. This can be used on any
        Windows machine without access to MSYS shells or compilers. See the
        READMEBinaryPackage.txt contained within the package for further details.


3. EXTERNAL PACKAGE REFERENCE
=============================

Name   : 7-Zip
Version: 9.20 (or above)
URL    : http://www.7-zip.org/a/7z920-x64.msi
For    : Extracting .zip and .7z files (can also be used to extract .tar.gz, .tar.bz2 etc.)

Name   : CMake
Version: 3.2.3 (or above)
URL    : http://www.cmake.org/files/v3.2/cmake-3.2.3-win32-x86.zip
For    : Required for building metis and parmetis third-party libraries.

Name   : Git
Version: 1.9.5
URL    : https://git-scm.com/download/win
For    : Version control system. Choose "Use Git from the Windows command prompt" and "Check out as-is. Commit as-is" in installer.

Name   : Java
Version: Version 8 Update 60
URL    : http://www.java.com
For    : Not essential

Name   : Microsoft Visual C++ 2010 Redistributable Package (x64)
Version: 2010
URL    : http://www.microsoft.com/en-us/download/details.aspx?id=14632
For    : msvcr100.dll required by OpenMPI

Name   : MinGW-w64
Version: 4.8.2 (or above)
URL    : http://sourceforge.net/projects/mingw-w64/files/Toolchains%20targetting%20Win64/Personal%20Builds/mingw-builds/4.8.2/threads-win32/seh/x86_64-4.8.2-release-win32-seh-rt_v3-rev3.7z/download
For    : Windows port of GCC C/C++ compiler (also includes Fortran compiler).

Name   : MSYS
Version: Latest
URL    : http://www.mingw.org/download/installer
For    : Provides Linux-like emulation shell and utilities (eg. find, grep, ls, make etc.) and tools (eg. bison, flex, m4, yacc etc.).
         Only the "mingw-developer-toolkit" and "msys-base" is required, DO NOT download ming32-* compiler packages (we use MinGW-w64 instead).
         Assuming you installed this to the default /c/MinGW directory, you will need to "cp -p make.exe gmake.exe" in /c/MinGW/msys/1.0/bin (unless have Strawberry Perl).

Name   : OpenMPI
Version: 1.6.1 (or above)
URL    : http://www.open-mpi.org/software/ompi/v1.6/downloads/OpenMPI_v1.6.1-1_win64.exe
For    : Bulding and running OpenMPI applications. Do not add to PATH in installer.

Name   : ParaView
Version: 4.3.1 (or above)
URL    : http://www.paraview.org/download/
For    : Graphically visualising foam results.

Name   : PExports
Version: 0.44 (or above)
URL    : http://sourceforge.net/projects/mingw/files/MinGW/Extension/pexports/pexports-0.46/pexports-0.46-mingw32-bin.tar.xz/download
For    : Extracting symbols from OpenMPI DLLs to pass to dlltool.exe (supplied in MinGW-w64 package). Move directory "bin" to "pexports-0.46" after unpacking.

Name   : Pyhton
Version: 2.7
URL    : https://www.python.org/download/releases/2.7/
For    : Not essential

Name   : Strawberry Perl
Version: 5.20.2.1 (or above)
URL    : http://strawberryperl.com/download/5.20.2.1/strawberry-perl-5.20.2.1-64bit.msi
For    : Running Perl scripts, "pexports.exe" and "gmake.exe" utilities

Name   : Subversion
Version: 1.8.13
URL    : http://sourceforge.net/projects/win32svn/files/1.8.13/
For    : ??????????

Name   : wget
Version: 1.11.4-1 (or above)
URL    : http://downloads.sourceforge.net/gnuwin32/wget-1.11.4-1-bin.zip, http://downloads.sourceforge.net/gnuwin32/wget-1.11.4-1-dep.zip, http://downloads.sourceforge.net/gnuwin32/wget-1.11.4-1-doc.zip
For    : Automatically downloading files from internet (eg. for use in automated build scripts)

git clone --branch nextRelease git://git.code.sf.net/p/openfoam-extend/foam-extend-3.1


4. EXAMPLE OUTPUT FROM SOURCING etc/bashrc.mingw
================================================

$ . etc/bashrc.mingw
Setting environment variables for user-defined installed system tools and utilities ...

Sourcing: /c/Users/UserName/Git/foam-extend-3.2/etc/bashrc
Sourcing: /c/Users/UserName/Git/foam-extend-3.2/etc/prefs.sh.mingw
Sourcing: /c/Users/UserName/Git/foam-extend-3.2/etc/settings.sh
    MESQUITE_DIR is initialized to: /c/Users/UserName/Git/foam-extend-3.2/ThirdParty/packages/mesquite-2.1.2
    METIS_DIR is initialized to: /c/Users/UserName/Git/foam-extend-3.2/ThirdParty/packages/metis-5.1.0
    PARMETIS_DIR is initialized to: /c/Users/UserName/Git/foam-extend-3.2/ThirdParty/packages/parmetis-4.0.3
    PARMGRIDGEN_DIR is initialized to: /c/Users/UserName/Git/foam-extend-3.2/ThirdParty/packages/ParMGridGen-1.0
    LIBCCMIO_DIR is initialized to:
    SCOTCH_DIR is initialized to: /c/Users/UserName/Git/foam-extend-3.2/ThirdParty/packages/scotch_6.0.0
    SCOTCH_DIR is initialized to: /c/Users/UserName/Git/foam-extend-3.2/ThirdParty/packages/scotch_6.0.0
    CMAKE_DIR is initialized to:
    M4_DIR is initialized to:
    BISON_DIR is initialized to:
    FLEX_DIR is initialized to: /C/MinGW/msys/1.0//bin/..
    ZOLTAN_DIR is initialized to:
    PYTHON_DIR is initialized to:
    PYFOAM_DIR is initialized to:
    PYFOAM_SITE_DIR is initialized to: /c/Users/UserName/Git/foam-extend-3.2/ThirdParty/PyFoamSiteScripts
    HWLOC_DIR is initialized to:
    QT_DIR is initialized to:
    PARAVIEW_DIR is initialized to:
    LLVM_DIR is initialized to:
    MESA_DIR is initialized to:
Sourcing: /c/Users/UserName/Git/foam-extend-3.2/etc/aliases.sh

Adding user-defined installed system tools to PATH ...
Setting OpenMPI environment settings ...

Checking versions of installed system tools (based on PATH) ...
7-Zip:       /c/Program Files/7-Zip/7z.exe [ 7-Zip [64] 9.20 Copyright (c) 1999-2010 Igor Pavlov 2010-11-18]
Bison:       /bin/bison.exe [bison (GNU Bison) 2.4.2]
CMake:       /c/Programs/cmake-3.2.3-win32-x86/bin/cmake.exe [cmake version 3.2.3]
Flex:        /bin/flex.exe [flex 2.5.35]
G++:         /c/Programs/mingw64/bin/g++.exe [g++.exe (x86_64-win32-seh-rev3, Built by MinGW-W64 project) 4.8.2]
GCC:         /c/Programs/mingw64/bin/gcc.exe [gcc.exe (x86_64-win32-seh-rev3, Built by MinGW-W64 project) 4.8.2]
GMake:       /c/Programs/strawberry-perl-5.20.2.1-64bit/c/bin/gmake.exe [GNU Make 4.0.90]
Git:         /c/Programs/Git/cmd/git.exe [git version 1.9.5.msysgit.1]
Java:        /c/ProgramData/Oracle/Java/javapath/java.exe [java version "1.8.0_45"]
M4:          /bin/m4.exe [m4 (GNU M4) 1.4.16]
Make:        /bin/make.exe [GNU Make 3.81]
MinGW-w64:   /c/Programs/mingw64
OpenMPI:     /c/Programs/OpenMPI_v1.6.1-x64/bin/mpirun.exe [mpirun.exe (OpenRTE) 1.6.1]
PEexports:   /c/Programs/pexports-0.46/bin/pexports.exe [PExports 0.44 Copyright 1998, Anders Norlander]
ParaView:    /c/Programs/ParaView-4.3.1
Perl:        /bin/perl.exe [ This is perl, v5.8.8 built for msys-64int]
Python:      /c/Programs/Python27/python.exe [Python 2.7.9]
Subversion:  /c/Programs/svn-win32-1.8.13/bin/svn.exe [svn, version 1.8.13 (r1667537)]
Vim:         /bin/vim.exe [VIM - Vi IMproved 7.3 (2010 Aug 15, compiled Mar 19 2011 15:37:04)]
Wget:        /c/Programs/wget-1.11.4-1/bin/wget.exe [GNU Wget 1.11.4 Copyright (C) 2008 Free Software Foundation, Inc.]


FOAM_INST_DIR=/c/Users/UserName/Git
WM_PROJECT_DIR=/c/Users/UserName/Git/foam-extend-3.2
WM_OSTYPE=MSWindows
ENVIRONMENT SETUP COMPLETE.
