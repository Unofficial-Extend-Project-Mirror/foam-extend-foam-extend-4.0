#!/bin/bash
#----------------------------------*-sh-*--------------------------------------
# =========                 |
# \\      /  F ield         | foam-extend: Open Source CFD
#  \\    /   O peration     | Version:     4.0
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
#     etc/settings.sh
#
# Description
#     Startup file for FOAM
#     Sourced from FOAM-??/etc/bashrc
#
#------------------------------------------------------------------------------

# prefix to PATH
_foamAddPath()
{
    while [ $# -ge 1 ]
    do
        export PATH=$1:$PATH
        shift
    done
}

# prefix to LD_LIBRARY_PATH
_foamAddLib()
{
    while [ $# -ge 1 ]
    do
        export LD_LIBRARY_PATH=$1:$LD_LIBRARY_PATH
        if [ "$WM_ARCH_BASE" = "darwin" ]
        then
            # do NOT add the lib of MacPort as this might break programs
            if [ "$1" != "/opt/local/lib" ]
            then
                export DYLD_LIBRARY_PATH=$1:$DYLD_LIBRARY_PATH
            fi
        fi
        shift
    done
}

# Source files, possibly with some verbosity
# Yes, this is the same definition as in the file etc/bash
# We need that definition available for scripts sourcing
# settings.sh directly.
_foamSource()
{
   while [ $# -ge 1 ]
   do
      [ "$FOAM_VERBOSE" -a "$PS1" ] && echo "Sourcing: $1"
      . $1
      shift
   done
}

# location of the jobControl directory
export FOAM_JOB_DIR=$HOME/$WM_PROJECT/jobControl

# wmake configuration
export WM_DIR=$WM_PROJECT_DIR/wmake
export WM_LINK_LANGUAGE=c++
export WM_OPTIONS=$WM_ARCH$WM_COMPILER$WM_PRECISION_OPTION$WM_COMPILE_OPTION
export PATH=$WM_DIR:$PATH

# base configuration
export FOAM_APP=$WM_PROJECT_DIR/applications
export FOAM_APPBIN=$WM_PROJECT_DIR/applications/bin/$WM_OPTIONS
export FOAM_LIB=$WM_PROJECT_DIR/lib
export FOAM_LIBBIN=$WM_PROJECT_DIR/lib/$WM_OPTIONS
export FOAM_SRC=$WM_PROJECT_DIR/src

# shared site configuration - similar naming convention as ~FOAM expansion
export FOAM_SITE_DIR=$WM_PROJECT_INST_DIR/site/$WM_PROJECT_VERSION
export FOAM_SITE_APPBIN=$FOAM_SITE_DIR/bin/$WM_OPTIONS
export FOAM_SITE_LIBBIN=$FOAM_SITE_DIR/lib/$WM_OPTIONS

# user configuration
export FOAM_USER_APPBIN=$WM_PROJECT_USER_DIR/applications/bin/$WM_OPTIONS
export FOAM_USER_LIBBIN=$WM_PROJECT_USER_DIR/lib/$WM_OPTIONS

# convenience
export FOAM_TUTORIALS=$WM_PROJECT_DIR/tutorials
export FOAM_UTILITIES=$FOAM_APP/utilities
export FOAM_SOLVERS=$FOAM_APP/solvers
export FOAM_RUN=$WM_PROJECT_USER_DIR/run
export FOAM_TEST_HARNESS_DIR=$WM_PROJECT_DIR/testHarness/foam-extend/$WM_PROJECT_VERSION

# add FOAM scripts and wmake to the path
export PATH=$WM_DIR:$WM_PROJECT_DIR/bin:$PATH

_foamAddPath $FOAM_APPBIN $FOAM_SITE_APPBIN $FOAM_USER_APPBIN
_foamAddLib  $FOAM_LIBBIN $FOAM_SITE_LIBBIN $FOAM_USER_LIBBIN


# Compiler settings
# ~~~~~~~~~~~~~~~~~
unset compilerBin compilerLib

# Select compiler installation
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# compilerInstall = FOAM | System
#
# We can override the value of compilerInstall from prefs.sh
: ${compilerInstall:=System}

# Or we can force it right here
#compilerInstall=FOAM
#compilerInstall=System

case "${compilerInstall}" in
FOAM)
    case "$WM_COMPILER" in
    Gcc)
        export WM_COMPILER_DIR=$WM_THIRD_PARTY_DIR/packages/gcc-4.6.4/platforms/$WM_OPTIONS
        _foamSource $WM_THIRD_PARTY_DIR/packages/gmp-5.1.2/platforms/$WM_OPTIONS/etc/gmp-5.1.2.sh
        _foamSource $WM_THIRD_PARTY_DIR/packages/mpfr-3.1.2/platforms/$WM_OPTIONS/etc/mpfr-3.1.2.sh
        _foamSource $WM_THIRD_PARTY_DIR/packages/mpc-1.0.1/platforms/$WM_OPTIONS/etc/mpc-1.0.1.sh
        _foamSource $WM_THIRD_PARTY_DIR/packages/gcc-4.6.4/platforms/$WM_OPTIONS/etc/gcc-4.6.4.sh
        ;;
    Gcc49)
        export WM_COMPILER_DIR=$WM_THIRD_PARTY_DIR/packages/gcc-4.9.2/platforms/$WM_OPTIONS
        _foamSource $WM_THIRD_PARTY_DIR/packages/gcc-4.9.2/platforms/$WM_OPTIONS/etc/gcc-4.9.2.sh
        ;;
    Gcc48)
        export WM_COMPILER_DIR=$WM_THIRD_PARTY_DIR/packages/gcc-4.8.4/platforms/$WM_OPTIONS
        _foamSource $WM_THIRD_PARTY_DIR/packages/gcc-4.8.4/platforms/$WM_OPTIONS/etc/gcc-4.8.4.sh
        ;;
    Gcc47)
        export WM_COMPILER_DIR=$WM_THIRD_PARTY_DIR/packages/gcc-4.7.4/platforms/$WM_OPTIONS
        _foamSource $WM_THIRD_PARTY_DIR/packages/gmp-5.1.2/platforms/$WM_OPTIONS/etc/gmp-5.1.2.sh
        _foamSource $WM_THIRD_PARTY_DIR/packages/mpfr-3.1.2/platforms/$WM_OPTIONS/etc/mpfr-3.1.2.sh
        _foamSource $WM_THIRD_PARTY_DIR/packages/mpc-1.0.1/platforms/$WM_OPTIONS/etc/mpc-1.0.1.sh
        _foamSource $WM_THIRD_PARTY_DIR/packages/gcc-4.7.4/platforms/$WM_OPTIONS/etc/gcc-4.7.4.sh
        ;;
    Gcc45)
        export WM_COMPILER_DIR=$WM_THIRD_PARTY_DIR/packages/gcc-4.5.1/platforms/$WM_OPTIONS
        _foamSource $WM_THIRD_PARTY_DIR/packages/gmp-5.0.1/platforms/$WM_OPTIONS/etc/gmp-5.0.1.sh
        _foamSource $WM_THIRD_PARTY_DIR/packages/mpfr-3.0.1/platforms/$WM_OPTIONS/etc/mpfr-3.0.1.sh
        _foamSource $WM_THIRD_PARTY_DIR/packages/mpc-0.8.2/platforms/$WM_OPTIONS/etc/mpc-0.8.2.sh
        _foamSource $WM_THIRD_PARTY_DIR/packages/gcc-4.5.1/platforms/$WM_OPTIONS/etc/gcc-4.5.1.sh
        ;;
    Gcc44)
        export WM_COMPILER_DIR=$WM_THIRD_PARTY_DIR/packages/gcc-4.4.5/platforms/$WM_OPTIONS
        _foamSource $WM_THIRD_PARTY_DIR/packages/gmp-5.0.1/platforms/$WM_OPTIONS/etc/gmp-5.0.1.sh
        _foamSource $WM_THIRD_PARTY_DIR/packages/mpfr-3.0.1/platforms/$WM_OPTIONS/etc/mpfr-3.0.1.sh
        _foamSource $WM_THIRD_PARTY_DIR/packages/gcc-4.4.5/platforms/$WM_OPTIONS/etc/gcc-4.4.5.sh
        ;;
    esac

    # Check that the compiler directory can be found
    if [ ! -d "$WM_COMPILER_DIR" ]
    then
        echo
        echo "Warning in $WM_PROJECT_DIR/etc/settings.sh:"
        echo "    Cannot find $WM_COMPILER_DIR installation."
        echo "    Please install this compiler version or if you wish to use the system compiler,"
        echo "    change the 'compilerInstall' setting to 'System' in this file"
        echo
    fi

    compilerBin=$WM_COMPILER_DIR/bin
    compilerLib=$WM_COMPILER_DIR/lib$WM_COMPILER_LIB_ARCH:$WM_COMPILER_DIR/lib
    ;;
esac

if [ -d "$compilerBin" ]
then
    _foamAddPath $compilerBin
    _foamAddLib  $compilerLib
fi

unset compilerBin compilerLib compilerInstall


if [ -z "$WM_CC" ]
then
    case "$WM_COMPILER" in
    Gcc*)
        export WM_CC='gcc'
        export WM_CXX='g++'
        ;;
    Icc)
        export WM_CC='icc'
        export WM_CXX='icpc'
        ;;
    esac
fi

# Communications library
# ~~~~~~~~~~~~~~~~~~~~~~

unset MPI_ARCH_PATH
mpi_version=unknown
case "$WM_MPLIB" in
OPENMPI)
    if [ ! -z $WM_THIRD_PARTY_USE_OPENMPI_188 ] && [ -e $WM_THIRD_PARTY_DIR/packages/openmpi-1.8.8/platforms/$WM_OPTIONS ]
        then
        mpi_version=openmpi-1.8.8
        if [ "$FOAM_VERBOSE" -a "$PS1" ]
        then
            echo "Using openmpi-1.8.8 from the ThirdParty package: $WM_THIRD_PARTY_DIR/packages/$mpi_version"
        fi
        _foamSource $WM_THIRD_PARTY_DIR/packages/$mpi_version/platforms/$WM_OPTIONS/etc/$mpi_version.sh

    elif [ ! -z $WM_THIRD_PARTY_USE_OPENMPI_165 ] && [ -e $WM_THIRD_PARTY_DIR/packages/openmpi-1.6.5/platforms/$WM_OPTIONS ]
        then
        mpi_version=openmpi-1.6.5
        if [ "$FOAM_VERBOSE" -a "$PS1" ]
        then
            echo "Using openmpi-1.6.5 from the ThirdParty package: $WM_THIRD_PARTY_DIR/packages/$mpi_version"
        fi
        _foamSource $WM_THIRD_PARTY_DIR/packages/$mpi_version/platforms/$WM_OPTIONS/etc/$mpi_version.sh

    elif [ ! -z $WM_THIRD_PARTY_USE_OPENMPI_143 ] && [ -e $WM_THIRD_PARTY_DIR/packages/openmpi-1.4.3/platforms/$WM_OPTIONS ]
        then
        mpi_version=openmpi-1.4.3
        if [ "$FOAM_VERBOSE" -a "$PS1" ]
        then
            echo "Using openmpi-1.4.3 from the ThirdParty package: $WM_THIRD_PARTY_DIR/packages/$mpi_version"
        fi
        _foamSource $WM_THIRD_PARTY_DIR/packages/$mpi_version/platforms/$WM_OPTIONS/etc/$mpi_version.sh

    elif [ ! -z $WM_THIRD_PARTY_USE_OPENMPI_15 ] && [ -e $WM_THIRD_PARTY_DIR/packages/openmpi-1.5/platforms/$WM_OPTIONS ]
        then
        mpi_version=openmpi-1.5
        if [ "$FOAM_VERBOSE" -a "$PS1" ]
        then
            echo "Using openmpi-1.5 from the ThirdParty package: $WM_THIRD_PARTY_DIR/packages/$mpi_version"
        fi
        _foamSource $WM_THIRD_PARTY_DIR/packages/$mpi_version/platforms/$WM_OPTIONS/etc/$mpi_version.sh
    fi

    # On Windows set mpi_version from value defined in bashrc.mingw:
    if [ "$WM_ARCH_BASE" == "mingw" ]
    then
        mpi_version=$MPI_VERSION_MINGW
    fi

    unset mpi_version
    ;;

MACPORTOPENMPI)
    unset OPAL_PREFIX

    export FOAM_MPI=openmpi-macport-$WM_MACPORT_MPI_VERSION

    # Currently not correctly working on MacPorts
    #    libDir=`mpicc-openmpi-$WM_MACPORT_MPI_VERSION --showme:libdirs`
    libDir=/opt/local/lib/openmpi-$WM_MACPORT_MPI_VERSION

    _foamAddLib     $libDir
    unset libDir

    which mpirun >/dev/null
    if [ $? -ne 0 ]
    then
        export WM_MPIRUN_PROG=mpirun-openmpi-$WM_MACPORT_MPI_VERSION
    fi
    ;;

MACPORTMPICH)
    export FOAM_MPI=mpich-macports-$WM_MACPORT_MPI_VERSION
    export MPI_HOME=$WM_THIRD_PARTY_DIR/$FOAM_MPI

    libDir=/opt/local/lib/mpich-$WM_MACPORT_MPI_VERSION

    _foamAddLib     $libDir
    unset libDir

    which mpirun >/dev/null
    if [ $? -ne 0 ]
    then
        export WM_MPIRUN_PROG=mpirun-mpich-$WM_MACPORT_MPI_VERSION
    fi
    ;;

SYSTEMOPENMPI)
    mpi_version=openmpi-system

    # make sure not the "old" mpi is used
    # Not sure if this is necessary anymore.
    # export OPAL_PREFIX=

    # Make sure OPENMPI_BIN_DIR is set and valid
    if [ -n "${OPENMPI_BIN_DIR}" ] && [ -d "${OPENMPI_BIN_DIR}" ]
    then
    # User defined value specified for OPENMPI_BIN_DIR
    #
    # WARNING:
    #          We assume this path specified by $OPENMPI_BIN_DIR is valid
    #          We assume the command mpicc is located somewhere under this path
    #          We assume the file mpi.h is located somewhere under this path
    #
    #          Otherwise, please double check your openmpi installation, you are
    #          probably missing the openmpi runtime and/or development packages
    #          available for your system.
    #
    _foamAddPath $OPENMPI_BIN_DIR
    else
    # Here, we assume your environment is already set for running
    # and developping with openmpi.
    #
    # Initialize OPENMPI_BIN_DIR using the path to mpicc
    export OPENMPI_BIN_DIR=$(dirname `which mpicc`)
    fi

    # Make sure OPENMPI_LIB_DIR is set
    if [ ! -n "${OPENMPI_LIB_DIR}" ]
    then
        # Initialize OPENMPI_LIB_DIR using the path to mpicc
        export OPENMPI_LIB_DIR="`mpicc --showme:libdirs`"
    fi

    # Make sure the dynamic libraries are accessible
    [   -n "${OPENMPI_LIB_DIR}" ]       && _foamAddLib $OPENMPI_LIB_DIR

    export MPI_HOME=`dirname $OPENMPI_BIN_DIR`
    export MPI_ARCH_PATH=$MPI_HOME
    export OPAL_PREFIX=$MPI_ARCH_PATH

    # We initialize the rest of the environment using mpicc --showme:
    [ ! -n "${OPENMPI_INCLUDE_DIR}" ]   && export OPENMPI_INCLUDE_DIR="`mpicc --showme:incdirs`"
    [ ! -n "${OPENMPI_COMPILE_FLAGS}" ] && export OPENMPI_COMPILE_FLAGS="`mpicc --showme:compile`"
    [ ! -n "${OPENMPI_LINK_FLAGS}" ]    && export OPENMPI_LINK_FLAGS="`mpicc --showme:link`"

    #
    # WARNING: We assume the file mpi.h will be available under the directories identified
    #          by the variable $OPENMPI_INCLUDE_DIR. Otherwise, please double check your
    #          system openmpi installation.

    # Set compilation flags here instead of in wmake/rules/../mplibSYSTEMOPENMPI
    export PINC="${OPENMPI_COMPILE_FLAGS}"
    export PLIBS="${OPENMPI_LINK_FLAGS}"

    # No longer needed, but we keep this as a reference, just in case...
    #libDir=`echo "$PLIBS" | sed -e 's/.*-L\([^ ]*\).*/\1/'`
    #_foamAddLib $libDir

    if [ "$FOAM_VERBOSE" -a "$PS1" ]
    then
        echo "  Environment variables defined for OpenMPI:"
        echo "    OPENMPI_BIN_DIR       : $OPENMPI_BIN_DIR"
        echo "    OPENMPI_LIB_DIR       : $OPENMPI_LIB_DIR"
        echo "    OPENMPI_INCLUDE_DIR   : $OPENMPI_INCLUDE_DIR"
        echo "    OPENMPI_COMPILE_FLAGS : $OPENMPI_COMPILE_FLAGS"
        echo "    OPENMPI_LINK_FLAGS    : $OPENMPI_LINK_FLAGS"
        echo ""
        echo "    MPI_HOME              : $MPI_HOME"
        echo "    MPI_ARCH_PATH         : $MPI_ARCH_PATH"
        echo "    OPAL_PREFIX           : $OPAL_PREFIX"
        echo "    PINC                  : $PINC"
        echo "    PLIBS                 : $PLIBS"
    fi

    unset mpi_version
    ;;


MVAPICH2)
    mpi_version=mvapich2-2.2
    if [ ! -z $WM_THIRD_PARTY_USE_MVAPICH2_20 ] && [ -e $WM_THIRD_PARTY_DIR/packages/mvapich2-2.2/platforms/$WM_OPTIONS ]
        then
        if [ "$FOAM_VERBOSE" -a "$PS1" ]
        then
            echo "Using mvapich2 from the ThirdParty package: $WM_THIRD_PARTY_DIR/packages/$mpi_version"
        fi
        _foamSource $WM_THIRD_PARTY_DIR/packages/$mpi_version/platforms/$WM_OPTIONS/etc/$mpi_version.sh
    fi

    export MPI_HOME=`dirname $MVAPICH2_BIN_DIR`
    export MPI_ARCH_PATH=$MPI_HOME

    export PINC="`mpicc -show -cc= -nativelinking`"
    export PLIBS="`mpicc -show -cc= | sed "s%$PINC%%"`"

    if [ "$FOAM_VERBOSE" -a "$PS1" ]
    then
        echo "  Environment variables defined for MVAPICH2:"
        echo "    MPI_ARCH_PATH         : $MPI_ARCH_PATH"
        echo "    PINC                  : $PINC"
        echo "    PLIBS                 : $PLIBS"
    fi

    unset mpi_version
    ;;

MPICH)
    mpi_version=mpich-1.2.4
    export MPI_HOME=$WM_THIRD_PARTY_DIR/$mpi_version
    export MPI_ARCH_PATH=$MPI_HOME/platforms/$WM_OPTIONS
    export MPICH_ROOT=$MPI_ARCH_PATH

    _foamAddPath $MPI_ARCH_PATH/bin
    _foamAddLib  $MPI_ARCH_PATH/lib

    unset mpi_version
    ;;

MPICH-GM)
    export MPI_ARCH_PATH=/opt/mpi
    export MPICH_PATH=$MPI_ARCH_PATH
    export MPICH_ROOT=$MPI_ARCH_PATH
    export GM_LIB_PATH=/opt/gm/lib64

    _foamAddPath $MPI_ARCH_PATH/bin
    _foamAddLib  $MPI_ARCH_PATH/lib
    _foamAddLib  $GM_LIB_PATH

    ;;

HPMPI)
    export MPI_HOME=/opt/hpmpi
    export MPI_ARCH_PATH=$MPI_HOME
    export MPICH_ROOT=$MPI_ARCH_PATH

    _foamAddPath $MPI_ARCH_PATH/bin

    case `uname -m` in
    i686)
        _foamAddLib $MPI_ARCH_PATH/lib/linux_ia32
        ;;

    x86_64)
        _foamAddLib $MPI_ARCH_PATH/lib/linux_amd64
        ;;
    ia64)
        _foamAddLib $MPI_ARCH_PATH/lib/linux_ia64
        ;;
    *)
        echo Unknown processor type `uname -m` for Linux
        ;;
    esac

    ;;

GAMMA)
    export MPI_ARCH_PATH=/usr
    ;;

MPI)
    export MPI_ARCH_PATH=/opt/mpi
    ;;

FJMPI)
    export MPI_ARCH_PATH=/opt/FJSVmpi2

    _foamAddPath $MPI_ARCH_PATH/bin
    _foamAddLib  $MPI_ARCH_PATH/lib/sparcv9
    _foamAddLib  /opt/FSUNf90/lib/sparcv9
    _foamAddLib  /opt/FJSVpnidt/lib
    ;;

QSMPI)
    export MPI_ARCH_PATH=/usr/lib/mpi

    _foamAddPath $MPI_ARCH_PATH/bin
    _foamAddLib $MPI_ARCH_PATH/lib
    ;;

SGIMPI)
    # no trailing slash
    [ "${MPI_ROOT%/}" = "${MPI_ROOT}" ] || MPI_ROOT="${MPI_ROOT%/}"

    export FOAM_MPI="${MPI_ROOT##*/}"
    export MPI_ARCH_PATH=$MPI_ROOT

    if [ ! -d "$MPI_ROOT" -o -z "$MPI_ARCH_PATH" ]
    then
        echo "Warning in $WM_PROJECT_DIR/etc/config/settings.sh:" 1>&2
        echo "    MPI_ROOT not a valid mpt installation directory or ending" \
             " in a '/'." 1>&2
        echo "    Please set MPI_ROOT to the mpt installation directory." 1>&2
        echo "    MPI_ROOT currently set to '$MPI_ROOT'" 1>&2
    fi

    if [ "$FOAM_VERBOSE" -a "$PS1" ]
    then
        echo "Using SGI MPT:" 1>&2
        echo "    MPI_ROOT : $MPI_ROOT" 1>&2
        echo "    FOAM_MPI : $FOAM_MPI" 1>&2
    fi

    _foamAddPath    $MPI_ARCH_PATH/bin
    _foamAddLib     $MPI_ARCH_PATH/lib

    export FOAM_MPI_LIBBIN=$FOAM_LIBBIN/sgimpi
    ;;

INTELMPI)
    # no trailing slash
    [ "${MPI_ROOT%/}" = "${MPI_ROOT}" ] || MPI_ROOT="${MPI_ROOT%/}"

    export FOAM_MPI="${MPI_ROOT##*/}"
    export MPI_ARCH_PATH=$MPI_ROOT

    if [ ! -d "$MPI_ROOT" -o -z "$MPI_ARCH_PATH" ]
    then
        echo "Warning in $WM_PROJECT_DIR/etc/config/settings.sh:" 1>&2
        echo "    MPI_ROOT not a valid mpt installation directory or ending" \
             " in a '/'." 1>&2
        echo "    Please set MPI_ROOT to the mpt installation directory." 1>&2
        echo "    MPI_ROOT currently set to '$MPI_ROOT'" 1>&2
    fi

    if [ "$FOAM_VERBOSE" -a "$PS1" ]
    then
        echo "Using INTEL MPI:" 1>&2
        echo "    MPI_ROOT : $MPI_ROOT" 1>&2
        echo "    FOAM_MPI : $FOAM_MPI" 1>&2
    fi

    _foamAddPath    $MPI_ARCH_PATH/intel64/bin
    _foamAddLib     $MPI_ARCH_PATH/intel64/lib

    export FOAM_MPI_LIBBIN=$FOAM_LIBBIN/intelmpi
    ;;

*)
    ;;
esac

# Removed $FOAM_MPI_LIBBIN.  HJ, 8/Aug/2015


# Set the minimum MPI buffer size (used by all platforms except SGI MPI)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
minBufferSize=20000000

if [ "${MPI_BUFFER_SIZE:=$minBufferSize}" -lt $minBufferSize ]
then
    MPI_BUFFER_SIZE=$minBufferSize
fi
export MPI_BUFFER_SIZE


# CUDA if available
# ~~~~~~~~~~~~~~~~~
[ -z "$CUDA_SYSTEM" ] && [ -e /usr/local/cuda-5.5/bin/nvcc ] && {
    export CUDA_DIR=/usr/local/cuda-5.5
    export CUDA_BIN_DIR=$CUDA_DIR/bin
    export CUDA_LIB_DIR=$CUDA_DIR/lib64
    export CUDA_INCLUDE_DIR=$CUDA_DIR/include
}

[ -d "$CUDA_LIB_DIR" ] && _foamAddPath $CUDA_BIN_DIR && _foamAddLib $CUDA_LIB_DIR


# CGAL library if available
# ~~~~~~~~~~~~~~~~~~~~~~~~~
[ -d "$CGAL_LIB_DIR" ] && _foamAddLib $CGAL_LIB_DIR


# Switch on the hoard memory allocator if available
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#if [ -f $FOAM_LIBBIN/libhoard.so ]
#then
#    export LD_PRELOAD=$FOAM_LIBBIN/libhoard.so:$LD_PRELOAD
#fi

# HOARD seems to crash Paraview on Mac OS X 10.6
#if [ -f $FOAM_LIBBIN/libhoard.dylib ]
#then
#    if [ -z "$DYLD_INSERT_LIBRARIES" ]
#    then
#        export DYLD_INSERT_LIBRARIES=$FOAM_LIBBIN/libhoard.dylib
#    else
#        export DYLD_INSERT_LIBRARIES=$FOAM_LIBBIN/libhoard.dylib:$DYLD_INSERT_LIBRARIES
#    fi
#fi


# Third party packages
#
# In order to use a pre-installed version of the ThirdParty packages, just set the
# appropriate XXX_SYSTEM environment variable for a given package in your prefs.sh
# file in order to disable the activation of the ThirdParty version of the same package.

# Load Mesquite library
# ~~~~~~~~~~~~~~~~~~~~~~
[ -z "$MESQUITE_SYSTEM" ] && [ ! -z $WM_THIRD_PARTY_USE_MESQUITE_212 ] && [ -e $WM_THIRD_PARTY_DIR/packages/mesquite-2.1.2/platforms/$WM_OPTIONS ] && {
    _foamSource $WM_THIRD_PARTY_DIR/packages/mesquite-2.1.2/platforms/$WM_OPTIONS/etc/mesquite-2.1.2.sh
}
[ "$FOAM_VERBOSE" -a "$PS1" ] && echo "    MESQUITE_DIR is initialized to: $MESQUITE_DIR"


# Load Metis library
# ~~~~~~~~~~~~~~~~~~
[ -z "$METIS_SYSTEM" ] && [ ! -z $WM_THIRD_PARTY_USE_METIS_510 ] && [ -e $WM_THIRD_PARTY_DIR/packages/metis-5.1.0/platforms/$WM_OPTIONS ] && {
    _foamSource $WM_THIRD_PARTY_DIR/packages/metis-5.1.0/platforms/$WM_OPTIONS/etc/metis-5.1.0.sh
}
[ "$FOAM_VERBOSE" -a "$PS1" ] && echo "    METIS_DIR is initialized to: $METIS_DIR"


# Load ParMetis library
# ~~~~~~~~~~~~~~~~~~~~~
[ -z "$PARMETIS_SYSTEM" ] && [ ! -z $WM_THIRD_PARTY_USE_PARMETIS_403 ] && [ -e $WM_THIRD_PARTY_DIR/packages/parmetis-4.0.3/platforms/$WM_OPTIONS ] && {
    _foamSource $WM_THIRD_PARTY_DIR/packages/parmetis-4.0.3/platforms/$WM_OPTIONS/etc/parmetis-4.0.3.sh
}
[ "$FOAM_VERBOSE" -a "$PS1" ] && echo "    PARMETIS_DIR is initialized to: $PARMETIS_DIR"


# Load ParMGridGen library
# ~~~~~~~~~~~~~~~~~~~~~~~~~
[ -z "$PARMGRIDGEN_SYSTEM" ] && [ ! -z $WM_THIRD_PARTY_USE_PARMGRIDGEN_10 ] && [ -e $WM_THIRD_PARTY_DIR/packages/ParMGridGen-1.0/platforms/$WM_OPTIONS ] && {
    _foamSource $WM_THIRD_PARTY_DIR/packages/ParMGridGen-1.0/platforms/$WM_OPTIONS/etc/ParMGridGen-1.0.sh
}
[ "$FOAM_VERBOSE" -a "$PS1" ] && echo "    PARMGRIDGEN_DIR is initialized to: $PARMGRIDGEN_DIR"


# Load Libccmio library
# ~~~~~~~~~~~~~~~~~~~~~
[ -z "$LIBCCMIO_SYSTEM" ] && [ ! -z $WM_THIRD_PARTY_USE_LIBCCMIO_261 ] && [ -e $WM_THIRD_PARTY_DIR/packages/libccmio-2.6.1/platforms/$WM_OPTIONS ] && {
    _foamSource $WM_THIRD_PARTY_DIR/packages/libccmio-2.6.1/platforms/$WM_OPTIONS/etc/libccmio-2.6.1.sh
}
[ "$FOAM_VERBOSE" -a "$PS1" ] && echo "    LIBCCMIO_DIR is initialized to: $LIBCCMIO_DIR"


# Load Scotch library
# ~~~~~~~~~~~~~~~~~~~
[ -z "$SCOTCH_SYSTEM" ] && [ ! -z $WM_THIRD_PARTY_USE_SCOTCH_604 ] && [ -e $WM_THIRD_PARTY_DIR/packages/scotch-6.0.4/platforms/$WM_OPTIONS ] && {
    _foamSource $WM_THIRD_PARTY_DIR/packages/scotch-6.0.4/platforms/$WM_OPTIONS/etc/scotch-6.0.4.sh
}
[ "$FOAM_VERBOSE" -a "$PS1" ] && echo "    SCOTCH_DIR is initialized to: $SCOTCH_DIR"

[ -z "$SCOTCH_SYSTEM" ] && [ ! -z $WM_THIRD_PARTY_USE_SCOTCH_600 ] && [ -e $WM_THIRD_PARTY_DIR/packages/scotch-6.0.0/platforms/$WM_OPTIONS ] && {
    _foamSource $WM_THIRD_PARTY_DIR/packages/scotch-6.0.0/platforms/$WM_OPTIONS/etc/scotch-6.0.0.sh
}
[ "$FOAM_VERBOSE" -a "$PS1" ] && echo "    SCOTCH_DIR is initialized to: $SCOTCH_DIR"


# Load cmake
# ~~~~~~~~~~
[ -z "$CMAKE_SYSTEM" ] && [ ! -z $WM_THIRD_PARTY_USE_CMAKE_322 ] && [ -e $WM_THIRD_PARTY_DIR/packages/cmake-3.2.2/platforms/$WM_OPTIONS ] && {
    _foamSource $WM_THIRD_PARTY_DIR/packages/cmake-3.2.2/platforms/$WM_OPTIONS/etc/cmake-3.2.2.sh
}
[ -z "$CMAKE_SYSTEM" ] && [ ! -z $WM_THIRD_PARTY_USE_CMAKE_2812 ] && [ -e $WM_THIRD_PARTY_DIR/packages/cmake-2.8.12/platforms/$WM_OPTIONS ] && {
    _foamSource $WM_THIRD_PARTY_DIR/packages/cmake-2.8.12/platforms/$WM_OPTIONS/etc/cmake-2.8.12.sh
}
[ "$FOAM_VERBOSE" -a "$PS1" ] && echo "    CMAKE_DIR is initialized to: $CMAKE_DIR"

# Load m4
# ~~~~~~~~~~
[ -z "$M4_SYSTEM" ] && [ ! -z $WM_THIRD_PARTY_USE_M4_1416 ] && [ -e $WM_THIRD_PARTY_DIR/packages/m4-1.4.16/platforms/$WM_OPTIONS ] && {
    _foamSource $WM_THIRD_PARTY_DIR/packages/m4-1.4.16/platforms/$WM_OPTIONS/etc/m4-1.4.16.sh
}
[ "$FOAM_VERBOSE" -a "$PS1" ] && echo "    M4_DIR is initialized to: $M4_DIR"

# Load bison
# ~~~~~~~~~~
[ -z "$BISON_SYSTEM" ] && [ ! -z $WM_THIRD_PARTY_USE_BISON_27 ] && [ -e $WM_THIRD_PARTY_DIR/packages/bison-2.7/platforms/$WM_OPTIONS ] && {
    _foamSource $WM_THIRD_PARTY_DIR/packages/bison-2.7/platforms/$WM_OPTIONS/etc/bison-2.7.sh
}
[ "$FOAM_VERBOSE" -a "$PS1" ] && echo "    BISON_DIR is initialized to: $BISON_DIR"

# Load flex
# ~~~~~~~~~~
[ -z "$FLEX_SYSTEM" ] && [ ! -z $WM_THIRD_PARTY_USE_FLEX_2535 ] && [ -e $WM_THIRD_PARTY_DIR/packages/flex-2.5.35/platforms/$WM_OPTIONS ] && {
    _foamSource $WM_THIRD_PARTY_DIR/packages/flex-2.5.35/platforms/$WM_OPTIONS/etc/flex-2.5.35.sh
}
[ "$FOAM_VERBOSE" -a "$PS1" ] && echo "    FLEX_DIR is initialized to: $FLEX_DIR"


# Load zoltan
# ~~~~~~~~~~
[ -z "$ZOLTAN_SYSTEM" ] && [ ! -z $WM_THIRD_PARTY_USE_ZOLTAN_36 ] && [ -e $WM_THIRD_PARTY_DIR/packages/zoltan_3.6 ] && {
    _foamSource $WM_THIRD_PARTY_DIR/packages/zoltan_3.6/platforms/$WM_OPTIONS/etc/zoltan_3.6.sh
}
[ "$FOAM_VERBOSE" -a "$PS1" ] && echo "    ZOLTAN_DIR is initialized to: $ZOLTAN_DIR"


# Load Python
# ~~~~~~~~~~~
[ -z "$PYTHON_SYSTEM" ] && [ ! -z $WM_THIRD_PARTY_USE_PYTHON_27 ] && [ -e $WM_THIRD_PARTY_DIR/packages/Python-2.7/platforms/$WM_OPTIONS ] && {
    _foamSource $WM_THIRD_PARTY_DIR/packages/Python-2.7/platforms/$WM_OPTIONS/etc/Python-2.7.sh
}
[ "$FOAM_VERBOSE" -a "$PS1" ] && echo "    PYTHON_DIR is initialized to: $PYTHON_DIR"

# Load PyFoam
# ~~~~~~~~~~~
[ -z "$PYFOAM_SYSTEM" ] && [ ! -z $WM_THIRD_PARTY_USE_PYFOAM_063 ] && [ -e $WM_THIRD_PARTY_DIR/packages/PyFoam-0.6.3/platforms/noarch ] && {
    _foamSource $WM_THIRD_PARTY_DIR/packages/PyFoam-0.6.3/platforms/noarch/etc/PyFoam-0.6.3.sh
}
[ -z "$PYFOAM_SYSTEM" ] && [ ! -z $WM_THIRD_PARTY_USE_PYFOAM_064 ] && [ -e $WM_THIRD_PARTY_DIR/packages/PyFoam-0.6.4/platforms/noarch ] && {
    _foamSource $WM_THIRD_PARTY_DIR/packages/PyFoam-0.6.4/platforms/noarch/etc/PyFoam-0.6.4.sh
}
[ "$FOAM_VERBOSE" -a "$PS1" ] && echo "    PYFOAM_DIR is initialized to: $PYFOAM_DIR"

# Make sure that there is a special directory for special PyFoam-scripts
[ -z "$PYFOAM_SITE_DIR" ] && export PYFOAM_SITE_DIR=$WM_THIRD_PARTY_DIR/PyFoamSiteScripts
if [ -d $PYFOAM_SITE_DIR/bin ]
then
    _foamAddPath $PYFOAM_SITE_DIR/bin
fi

[ "$FOAM_VERBOSE" -a "$PS1" ] && echo "    PYFOAM_SITE_DIR is initialized to: $PYFOAM_SITE_DIR"

# Load hwloc
# ~~~~~~~~~~~
[ -z "$HWLOC_SYSTEM" ] && [ ! -z $WM_THIRD_PARTY_USE_HWLOC_1101 ] && [ -e $WM_THIRD_PARTY_DIR/packages/hwloc-1.10.1/platforms/$WM_OPTIONS ] && {
    _foamSource $WM_THIRD_PARTY_DIR/packages/hwloc-1.10.1/platforms/$WM_OPTIONS/etc/hwloc-1.10.1.sh
}
[ -z "$HWLOC_SYSTEM" ] && [ ! -z $WM_THIRD_PARTY_USE_HWLOC_172 ] && [ -e $WM_THIRD_PARTY_DIR/packages/hwloc-1.7.2/platforms/$WM_OPTIONS ] && {
    _foamSource $WM_THIRD_PARTY_DIR/packages/hwloc-1.7.2/platforms/$WM_OPTIONS/etc/hwloc-1.7.2.sh
}
[ "$FOAM_VERBOSE" -a "$PS1" ] && echo "    HWLOC_DIR is initialized to: $HWLOC_DIR"

# Load QT
# ~~~~~~~
[ ! -z "$QT_THIRD_PARTY" ] && [ ! -z $WM_THIRD_PARTY_USE_QT_486 ] && [ -e $WM_THIRD_PARTY_DIR/packages/qt-everywhere-opensource-src-4.8.6/platforms/$WM_OPTIONS ] && {
    _foamSource $WM_THIRD_PARTY_DIR/packages/qt-everywhere-opensource-src-4.8.6/platforms/$WM_OPTIONS/etc/qt-everywhere-opensource-src-4.8.6.sh
}
[ "$FOAM_VERBOSE" -a "$PS1" ] && echo "    QT_DIR is initialized to: $QT_DIR"


# Load ParaView
# ~~~~~~~~~~~~~
[ -z "$PARAVIEW_SYSTEM" ] && [ ! -z $WM_THIRD_PARTY_USE_PARAVIEW_440 ] && [ -e $WM_THIRD_PARTY_DIR/packages/ParaView-4.4.0/platforms/$WM_OPTIONS ] && {
    _foamSource $WM_THIRD_PARTY_DIR/packages/ParaView-4.4.0/platforms/$WM_OPTIONS/etc/ParaView-4.4.0.sh
}
[ -z "$PARAVIEW_SYSTEM" ] && [ ! -z $WM_THIRD_PARTY_USE_PARAVIEW_431 ] && [ -e $WM_THIRD_PARTY_DIR/packages/ParaView-4.3.1/platforms/$WM_OPTIONS ] && {
    _foamSource $WM_THIRD_PARTY_DIR/packages/ParaView-4.3.1/platforms/$WM_OPTIONS/etc/ParaView-4.3.1.sh
}
[ -z "$PARAVIEW_SYSTEM" ] && [ ! -z $WM_THIRD_PARTY_USE_PARAVIEW_410 ] && [ -e $WM_THIRD_PARTY_DIR/packages/ParaView-4.1.0/platforms/$WM_OPTIONS ] && {
    _foamSource $WM_THIRD_PARTY_DIR/packages/ParaView-4.1.0/platforms/$WM_OPTIONS/etc/ParaView-4.1.0.sh
}
[ -z "$PARAVIEW_SYSTEM" ] && [ ! -z $WM_THIRD_PARTY_USE_PARAVIEW_401 ] && [ -e $WM_THIRD_PARTY_DIR/packages/ParaView-4.0.1/platforms/$WM_OPTIONS ] && {
    _foamSource $WM_THIRD_PARTY_DIR/packages/ParaView-4.0.1/platforms/$WM_OPTIONS/etc/ParaView-4.0.1.sh
}
[ -z "$PARAVIEW_SYSTEM" ] && [ ! -z $WM_THIRD_PARTY_USE_PARAVIEW_3141 ] && [ -e $WM_THIRD_PARTY_DIR/packages/ParaView-3.14.1/platforms/$WM_OPTIONS ] && {
    _foamSource $WM_THIRD_PARTY_DIR/packages/ParaView-3.14.1/platforms/$WM_OPTIONS/etc/ParaView-3.14.1.sh
}
[ -z "$PARAVIEW_SYSTEM" ] && [ ! -z $WM_THIRD_PARTY_USE_PARAVIEW_381 ] && [ -e $WM_THIRD_PARTY_DIR/packages/ParaView-3.8.1/platforms/$WM_OPTIONS ] && {
    _foamSource $WM_THIRD_PARTY_DIR/packages/ParaView-3.8.1/platforms/$WM_OPTIONS/etc/ParaView-3.8.1.sh
}
[ "$FOAM_VERBOSE" -a "$PS1" ] && echo "    PARAVIEW_DIR is initialized to: $PARAVIEW_DIR"


# Load llvm
# ~~~~~~~~~~~~~~~~~~~~~~
[ -z "$LLVM_SYSTEM" ] && [ ! -z $WM_THIRD_PARTY_USE_LLVM_360 ] && [ -e $WM_THIRD_PARTY_DIR/packages/llvm-3.6.0/platforms/$WM_OPTIONS ] && {
    _foamSource $WM_THIRD_PARTY_DIR/packages/llvm-3.6.0/platforms/$WM_OPTIONS/etc/llvm-3.6.0.sh
}
[ "$FOAM_VERBOSE" -a "$PS1" ] && echo "    LLVM_DIR is initialized to: $LLVM_DIR"

# Load mesa
# ~~~~~~~~~~~~~~~~~~~~~~
[ -z "$MESA_SYSTEM" ] && [ ! -z $WM_THIRD_PARTY_USE_MESA_GIT ] && [ -e $WM_THIRD_PARTY_DIR/packages/mesa-git/platforms/$WM_OPTIONS ] && {
    _foamSource $WM_THIRD_PARTY_DIR/packages/mesa-git/platforms/$WM_OPTIONS/etc/mesa-git.sh
}
[ "$FOAM_VERBOSE" -a "$PS1" ] && echo "    MESA_DIR is initialized to: $MESA_DIR"


# cleanup environment:
# ~~~~~~~~~~~~~~~~~~~~
unset _foamAddPath _foamAddLib minBufferSize

# -----------------------------------------------------------------------------
