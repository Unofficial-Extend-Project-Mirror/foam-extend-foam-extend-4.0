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
#     etc/settings.csh
#
# Description
#     Startup file for FOAM
#     Sourced from FOAM-??/etc/cshrc
#
#------------------------------------------------------------------------------

# prefix to PATH
alias _foamAddPath 'set path=(\!* $path)'
# prefix to LD_LIBRARY_PATH
alias _foamAddLib 'setenv LD_LIBRARY_PATH \!*\:${LD_LIBRARY_PATH}'

# location of the jobControl directory
setenv FOAM_JOB_DIR $WM_PROJECT_INST_DIR/jobControl

# wmake configuration
setenv WM_DIR $WM_PROJECT_DIR/wmake
setenv WM_LINK_LANGUAGE c++
setenv WM_OPTIONS $WM_ARCH$WM_COMPILER$WM_PRECISION_OPTION$WM_COMPILE_OPTION
set path=($WM_DIR $path)

# base configuration
setenv FOAM_APP $WM_PROJECT_DIR/applications
setenv FOAM_APPBIN $WM_PROJECT_DIR/applications/bin/$WM_OPTIONS
setenv FOAM_LIB $WM_PROJECT_DIR/lib
setenv FOAM_LIBBIN $WM_PROJECT_DIR/lib/$WM_OPTIONS
setenv FOAM_SRC $WM_PROJECT_DIR/src

# shared site configuration - similar naming convention as ~FOAM expansion
setenv FOAM_SITE_DIR $WM_PROJECT_INST_DIR/site/$WM_PROJECT_VERSION
setenv FOAM_SITE_APPBIN $FOAM_SITE_DIR/bin/$WM_OPTIONS
setenv FOAM_SITE_LIBBIN $FOAM_SITE_DIR/lib/$WM_OPTIONS

# user configuration
setenv FOAM_USER_APPBIN $WM_PROJECT_USER_DIR/applications/bin/$WM_OPTIONS
setenv FOAM_USER_LIBBIN $WM_PROJECT_USER_DIR/lib/$WM_OPTIONS

# convenience
setenv FOAM_TUTORIALS $WM_PROJECT_DIR/tutorials
setenv FOAM_UTILITIES $FOAM_APP/utilities
setenv FOAM_SOLVERS $FOAM_APP/solvers
setenv FOAM_RUN $WM_PROJECT_USER_DIR/run
setenv FOAM_TEST_HARNESS_DIR $WM_PROJECT_DIR/testHarness/foam-extend/$WM_PROJECT_VERSION

# add FOAM scripts and wmake to the path
set path=($WM_DIR $WM_PROJECT_DIR/bin $path)

_foamAddPath $FOAM_APPBIN
_foamAddPath $FOAM_SITE_APPBIN
_foamAddPath $FOAM_USER_APPBIN
 # Make sure to pick up dummy versions of external libraries last
_foamAddLib  $FOAM_LIBBIN/dummy
_foamAddLib  $FOAM_LIBBIN
_foamAddLib  $FOAM_SITE_LIBBIN
_foamAddLib  $FOAM_USER_LIBBIN


# Select compiler installation
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# compilerInstall = FOAM | System
#set compilerInstall=FOAM
#set compilerInstall=System
if ( ! $?compilerInstall ) then
    setenv compilerInstall System
endif

switch ("$compilerInstall")
case FOAM:
    switch ("$WM_COMPILER")
    case Gcc:
        setenv WM_COMPILER_DIR $WM_THIRD_PARTY_DIR/packages/gcc-4.6.4/platforms/$WM_OPTIONS
        _foamSource $WM_THIRD_PARTY_DIR/packages/gmp-5.1.2/platforms/$WM_OPTIONS/etc/gmp-5.1.2.csh
        _foamSource $WM_THIRD_PARTY_DIR/packages/mpfr-3.1.2/platforms/$WM_OPTIONS/etc/mpfr-3.1.2.csh
        _foamSource $WM_THIRD_PARTY_DIR/packages/mpc-1.0.1/platforms/$WM_OPTIONS/etc/mpc-1.0.1.csh
        _foamSource $WM_THIRD_PARTY_DIR/packages/gcc-4.6.4/platforms/$WM_OPTIONS/etc/gcc-4.6.4.csh
    breaksw
    case Gcc49:
        setenv WM_COMPILER_DIR $WM_THIRD_PARTY_DIR/packages/gcc-4.9.2/platforms/$WM_OPTIONS
        _foamSource $WM_THIRD_PARTY_DIR/packages/gcc-4.9.2/platforms/$WM_OPTIONS/etc/gcc-4.9.2.csh
    breaksw
    case Gcc48:
        setenv WM_COMPILER_DIR $WM_THIRD_PARTY_DIR/packages/gcc-4.8.4/platforms/$WM_OPTIONS
        _foamSource $WM_THIRD_PARTY_DIR/packages/gcc-4.8.4/platforms/$WM_OPTIONS/etc/gcc-4.8.4.csh
    breaksw
    case Gcc47:
        setenv WM_COMPILER_DIR $WM_THIRD_PARTY_DIR/packages/gcc-4.7.4/platforms/$WM_OPTIONS
        _foamSource $WM_THIRD_PARTY_DIR/packages/gmp-5.1.2/platforms/$WM_OPTIONS/etc/gmp-5.1.2.csh
        _foamSource $WM_THIRD_PARTY_DIR/packages/mpfr-3.1.2/platforms/$WM_OPTIONS/etc/mpfr-3.1.2.csh
        _foamSource $WM_THIRD_PARTY_DIR/packages/mpc-1.0.1/platforms/$WM_OPTIONS/etc/mpc-1.0.1.csh
        _foamSource $WM_THIRD_PARTY_DIR/packages/gcc-4.7.4/platforms/$WM_OPTIONS/etc/gcc-4.7.4.csh
    breaksw
    case Gcc45:
        setenv WM_COMPILER_DIR $WM_THIRD_PARTY_DIR/packages/gcc-4.5.1/platforms/$WM_OPTIONS
        _foamSource $WM_THIRD_PARTY_DIR/packages/gmp-5.0.1/platforms/$WM_OPTIONS/etc/gmp-5.0.1.csh
        _foamSource $WM_THIRD_PARTY_DIR/packages/mpfr-3.0.1/platforms/$WM_OPTIONS/etc/mpfr-3.0.1.csh
        _foamSource $WM_THIRD_PARTY_DIR/packages/mpc-0.8.2/platforms/$WM_OPTIONS/etc/mpc-0.8.2.csh
        _foamSource $WM_THIRD_PARTY_DIR/packages/gcc-4.5.1/platforms/$WM_OPTIONS/etc/gcc-4.5.1.csh
    breaksw
    case Gcc44:
        setenv WM_COMPILER_DIR $WM_THIRD_PARTY_DIR/packages/gcc-4.4.5/platforms/$WM_OPTIONS
        _foamSource $WM_THIRD_PARTY_DIR/packages/gmp-5.0.1/platforms/$WM_OPTIONS/etc/gmp-5.0.1.csh
        _foamSource $WM_THIRD_PARTY_DIR/packages/mpfr-3.0.1/platforms/$WM_OPTIONS/etc/mpfr-3.0.1.csh
        _foamSource $WM_THIRD_PARTY_DIR/packages/gcc-4.4.5/platforms/$WM_OPTIONS/etc/gcc-4.4.5.csh
    breaksw
    endsw

    # Check that the compiler directory can be found
    if ( ! -d "$WM_COMPILER_DIR" ) then
        echo
        echo "Warning in $WM_PROJECT_DIR/etc/settings.csh:"
        echo "    Cannot find $WM_COMPILER_DIR installation."
        echo "    Please install this compiler version or if you wish to use the system compiler,"
        echo "    change the 'compilerInstall' setting to 'System' in this file"
        echo
    endif

    _foamAddPath ${WM_COMPILER_DIR}/bin
    _foamAddLib  ${WM_COMPILER_DIR}/lib${WM_COMPILER_LIB_ARCH}
    _foamAddLib  ${WM_COMPILER_DIR}/lib

    breaksw
endsw


switch ("$WM_COMPILER")
case Gcc*:
    setenv WM_CC 'gcc'
    setenv WM_CXX 'g++'
    breaksw
case Icc:
    setenv WM_CC 'icc'
    setenv WM_CXX 'icpc'
    breaksw
endsw

# Communications library
# ~~~~~~~~~~~~~~~~~~~~~~

unset MPI_ARCH_PATH

set mpi_version=unknown

switch ("$WM_MPLIB")
case OPENMPI:
    if ($?WM_THIRD_PARTY_USE_OPENMPI_188 != 0 && -d $WM_THIRD_PARTY_DIR/packages/openmpi-1.8.8/platforms/$WM_OPTIONS ) then
        set mpi_version=openmpi-1.8.8

        if ($?FOAM_VERBOSE && $?prompt) then
            echo "Using openmpi-1.8.8 from the ThirdParty package: $WM_THIRD_PARTY_DIR/packages/$mpi_version"
        endif
        _foamSource $WM_THIRD_PARTY_DIR/packages/$mpi_version/platforms/$WM_OPTIONS/etc/$mpi_version.csh

    else if ($?WM_THIRD_PARTY_USE_OPENMPI_165 != 0 && -d $WM_THIRD_PARTY_DIR/packages/openmpi-1.6.5/platforms/$WM_OPTIONS ) then
        set mpi_version=openmpi-1.6.5

        if ($?FOAM_VERBOSE && $?prompt) then
            echo "Using openmpi-1.6.5 from the ThirdParty package: $WM_THIRD_PARTY_DIR/packages/$mpi_version"
        endif
        _foamSource $WM_THIRD_PARTY_DIR/packages/$mpi_version/platforms/$WM_OPTIONS/etc/$mpi_version.csh

    else if ($?WM_THIRD_PARTY_USE_OPENMPI_143 != 0 && -d $WM_THIRD_PARTY_DIR/packages/openmpi-1.4.3/platforms/$WM_OPTIONS ) then
        set mpi_version=openmpi-1.4.3

        if ($?FOAM_VERBOSE && $?prompt) then
            echo "Using openmpi-1.4.3 from the ThirdParty package: $WM_THIRD_PARTY_DIR/packages/$mpi_version"
        endif
        _foamSource $WM_THIRD_PARTY_DIR/packages/$mpi_version/platforms/$WM_OPTIONS/etc/$mpi_version.csh

    else if ($?WM_THIRD_PARTY_USE_OPENMPI_143 != 0 && -d $WM_THIRD_PARTY_DIR/packages/openmpi-1.5/platforms/$WM_OPTIONS ) then
        set mpi_version=openmpi-1.5

        if ($?FOAM_VERBOSE && $?prompt) then
            echo "Using openmpi-1.5 from the ThirdParty package: $WM_THIRD_PARTY_DIR/packages/$mpi_version"
        endif
        _foamSource $WM_THIRD_PARTY_DIR/packages/$mpi_version/platforms/$WM_OPTIONS/etc/$mpi_version.csh
    endif

    unset mpi_version
    breaksw

case SYSTEMOPENMPI:
    set mpi_version=openmpi-system

    # make sure not the "old" mpi is used
    # Not sure if this is necessary anymore.
    # export OPAL_PREFIX=

    # Make sure OPENMPI_BIN_DIR is set and valid
    if ($?OPENMPI_BIN_DIR != 0 ) then
        if (-d "${OPENMPI_BIN_DIR}" ) then
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
         endif
    else
    # Here, we assume your environment is already set for running
    # and developing with openmpi.
    #
    # Initialize OPENMPI_BIN_DIR using the path to mpicc
        set mpicc_cmd=`which mpicc`
        setenv OPENMPI_BIN_DIR `dirname $mpicc_cmd`
        unset mpicc_cmd
    endif

    # Make sure OPENMPI_LIB_DIR is set
    if ( $?OPENMPI_LIB_DIR == 0 ) then
        # Initialize OPENMPI_LIB_DIR using the path to mpicc
        setenv OPENMPI_LIB_DIR "`mpicc --showme:libdirs`"
    endif

    # Make sure the dynamic libraries are accessible
    if( $?OPENMPI_LIB_DIR != 0 ) then
        _foamAddLib $OPENMPI_LIB_DIR
    endif

    setenv MPI_HOME `dirname $OPENMPI_BIN_DIR`
    setenv MPI_ARCH_PATH $MPI_HOME
    setenv OPAL_PREFIX $MPI_ARCH_PATH

    # We initialize the rest of the environment using mpicc --showme:
    if ($?OPENMPI_INCLUDE_DIR == 0) then
        setenv OPENMPI_INCLUDE_DIR "`mpicc --showme:incdirs`"
    endif

    if (${?OPENMPI_COMPILE_FLAGS} == 0) then
        setenv OPENMPI_COMPILE_FLAGS "`mpicc --showme:compile`"
    endif

    if (${?OPENMPI_COMPILE_FLAGS} == 0) then
        setenv OPENMPI_COMPILE_FLAGS "`mpicc --showme:compile`"
    endif

    if (${?OPENMPI_LINK_FLAGS} == 0) then
        setenv OPENMPI_LINK_FLAGS "`mpicc --showme:link`"
    endif

    #
    # WARNING: We assume the file mpi.h will be available under the directories identified
    #          by the variable $OPENMPI_INCLUDE_DIR. Otherwise, please double check your
    #          system openmpi installation.

    # Set compilation flags here instead of in wmake/rules/../mplibSYSTEMOPENMPI
    setenv PINC "$OPENMPI_COMPILE_FLAGS"
    setenv PLIBS "$OPENMPI_LINK_FLAGS"

    # No longer needed, but we keep this as a reference, just in case...
    #libDir `echo "$PLIBS" | sed -e 's/.*-L\([^ ]*\).*/\1/'`
    #_foamAddLib $libDir

    if ($?FOAM_VERBOSE && $?prompt) then
        echo "Using system installed OpenMPI:"
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
    endif

    unset mpi_version
    breaksw

case MVAPICH2:
    set mpi_version=mvapich2

    if ($?MVAPICH2_BIN_DIR != 0) then
        if (-d "${MVAPICH2_BIN_DIR}" ) then
        _foamAddPath $MVAPICH2_BIN_DIR
        endif
    else
        set mpicc_cmd=`which mpicc`
        setenv MVAPICH2_BIN_DIR `dirname $mpicc_cmd`
        unset mpicc_cmd
    endif

    setenv MPI_HOME `dirname $MVAPICH2_BIN_DIR`
    setenv MPI_ARCH_PATH $MPI_HOME

    setenv  PINC "`mpicc -show -cc= -nativelinking`"
    setenv  PLIBS "`mpicc -show -cc= | sed "s%$PINC%%"`"

    if ($?FOAM_VERBOSE && $?prompt) then
        echo "  Environment variables defined for MVAPICH2:"
        echo "    MPI_ARCH_PATH         : $MPI_ARCH_PATH"
        echo "    PINC                  : $PINC"
        echo "    PLIBS                 : $PLIBS"
    endif

    unset mpi_version
    breaksw

case MPICH:
    set mpi_version=mpich-1.2.4
    setenv MPI_HOME $WM_THIRD_PARTY_DIR/$mpi_version
    setenv MPI_ARCH_PATH $MPI_HOME/platforms/$WM_OPTIONS
    setenv MPICH_ROOT $MPI_ARCH_PATH

    _foamAddPath $MPI_ARCH_PATH/bin
    _foamAddLib  $MPI_ARCH_PATH/lib

    unset mpi_version
    breaksw

case MPICH-GM:
    setenv MPI_ARCH_PATH /opt/mpi
    setenv MPICH_PATH $MPI_ARCH_PATH
    setenv MPICH_ROOT $MPI_ARCH_PATH
    setenv GM_LIB_PATH /opt/gm/lib64

    _foamAddPath $MPI_ARCH_PATH/bin
    _foamAddLib  $MPI_ARCH_PATH/lib
    _foamAddLib  $GM_LIB_PATH
    breaksw

case HPMPI:
    setenv MPI_HOME /opt/hpmpi
    setenv MPI_ARCH_PATH $MPI_HOME
    setenv MPICH_ROOT=$MPI_ARCH_PATH

    _foamAddPath $MPI_ARCH_PATH/bin

    switch (`uname -m`)
    case i686:
        _foamAddLib $MPI_ARCH_PATH/lib/linux_ia32
        breaksw
    case x86_64:
        _foamAddLib $MPI_ARCH_PATH/lib/linux_amd64
        breaksw
    case ia64:
        _foamAddLib $MPI_ARCH_PATH/lib/linux_ia64
        breaksw
    default:
        echo Unknown processor type `uname -m` for Linux
        breaksw
    endsw
    breaksw

case GAMMA:
    setenv MPI_ARCH_PATH /usr
    breaksw

case MPI:
    setenv MPI_ARCH_PATH /opt/mpi
    breaksw

case FJMPI:
    setenv MPI_ARCH_PATH /opt/FJSVmpi2
    _foamAddPath $MPI_ARCH_PATH/bin
    _foamAddLib  $MPI_ARCH_PATH/lib/sparcv9
    _foamAddLib  /opt/FSUNf90/lib/sparcv9
    _foamAddLib  /opt/FJSVpnidt/lib
    breaksw

case QSMPI:
    setenv MPI_ARCH_PATH /usr/lib/mpi
    _foamAddPath $MPI_ARCH_PATH/bin
    _foamAddLib $MPI_ARCH_PATH/lib

    breaksw

case SGIMPI:
    if ( ! $?MPI_ROOT) setenv MPI_ROOT /dummy

    if ( ! -d "$MPI_ROOT" ) then
        echo "Warning in $WM_PROJECT_DIR/etc/config/settings.csh:"
        echo "    MPI_ROOT not a valid mpt installation directory."
        echo "    Please set MPI_ROOT to the mpt installation directory."
        echo "    (usually done by loading the mpt module)"
        echo "    MPI_ROOT currently set to '$MPI_ROOT'"
    endif

    if ( "${MPI_ROOT:h}/" == $MPI_ROOT ) then
        setenv MPI_ROOT ${MPI_ROOT:h}
    endif

    setenv FOAM_MPI ${MPI_ROOT:t}
    setenv MPI_ARCH_PATH $MPI_ROOT


    if ($?FOAM_VERBOSE && $?prompt) then
        echo "Using SGI MPT:"
        echo "    MPI_ROOT : $MPI_ROOT"
        echo "    FOAM_MPI : $FOAM_MPI"
    endif


    _foamAddPath    $MPI_ARCH_PATH/bin
    _foamAddLib     $MPI_ARCH_PATH/lib
    breaksw

case INTELMPI:
    if ( ! $?MPI_ROOT) setenv MPI_ROOT /dummy

    if ( ! -d "$MPI_ROOT" ) then
        echo "Warning in $WM_PROJECT_DIR/etc/config/settings.csh:"
        echo "    MPI_ROOT not a valid mpt installation directory."
        echo "    Please set MPI_ROOT to the mpt installation directory."
        echo "    (usually done by loading the mpt module)"
        echo "    MPI_ROOT currently set to '$MPI_ROOT'"
    endif

    if ( "${MPI_ROOT:h}/" == $MPI_ROOT ) then
        setenv MPI_ROOT ${MPI_ROOT:h}
    endif

    setenv FOAM_MPI ${MPI_ROOT:t}
    setenv MPI_ARCH_PATH $MPI_ROOT


    if ($?FOAM_VERBOSE && $?prompt) then
        echo "Using INTEL MPT:"
        echo "    MPI_ROOT : $MPI_ROOT"
        echo "    FOAM_MPI : $FOAM_MPI"
    endif


    _foamAddPath    $MPI_ARCH_PATH/bin64
    _foamAddLib     $MPI_ARCH_PATH/lib64
    breaksw

default:
    breaksw
endsw

# Set the minimum MPI buffer size (used by all platforms except SGI MPI)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
set minBufferSize=20000000

if ( $?MPI_BUFFER_SIZE ) then
    if ( $MPI_BUFFER_SIZE < $minBufferSize ) then
        setenv MPI_BUFFER_SIZE $minBufferSize
    endif
else
    setenv MPI_BUFFER_SIZE $minBufferSize
endif

# CUDA library
# ~~~~~~~~~~~~
if ( $?CUDA_SYSTEM == 0 && -e /usr/local/cuda-5.5/bin/nvcc ) then
    setenv CUDA_DIR /usr/local/cuda-5.5
    setenv CUDA_BIN_DIR $CUDA_DIR/bin
    setenv CUDA_LIB_DIR $CUDA_DIR/lib64
    setenv CUDA_INCLUDE_DIR $CUDA_DIR/include
    setenv CUDA_ARCH sm_20
endif

if ( $?CUDA_BIN_DIR ) then
    _foamAddPath $CUDA_BIN_DIR
endif

if ( $?CUDA_LIB_DIR ) then
    _foamAddLib $CUDA_LIB_DIR
endif

# CGAL library if available
# ~~~~~~~~~~~~~~~~~~~~~~~~~
if ( $?CGAL_LIB_DIR ) then
    _foamAddLib $CGAL_LIB_DIR
endif

# Mesquite library if available
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if ( $?MESQUITE_SYSTEM == 0 && $?WM_THIRD_PARTY_USE_MESQUITE_212 != 0 && -e $WM_THIRD_PARTY_DIR/packages/mesquite-2.1.2/platforms/$WM_OPTIONS ) then
    _foamSource $WM_THIRD_PARTY_DIR/packages/mesquite-2.1.2/platforms/$WM_OPTIONS/etc/mesquite-2.1.2.csh
endif

# Metis library if available
# ~~~~~~~~~~~~~~~~~~~~~~~~~~
if ( $?METIS_SYSTEM == 0 && $?WM_THIRD_PARTY_USE_METIS_510 != 0 && -e $WM_THIRD_PARTY_DIR/packages/metis-5.1.0/platforms/$WM_OPTIONS ) then
    _foamSource $WM_THIRD_PARTY_DIR/packages/metis-5.1.0/platforms/$WM_OPTIONS/etc/metis-5.1.0.csh
endif

# ParMetis library if available
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if ( $?PARMETIS_SYSTEM == 0 && $?WM_THIRD_PARTY_USE_PARMETIS_403 != 0 && -e $WM_THIRD_PARTY_DIR/packages/parmetis-4.0.3/platforms/$WM_OPTIONS ) then
    _foamSource $WM_THIRD_PARTY_DIR/packages/parmetis-4.0.3/platforms/$WM_OPTIONS/etc/parmetis-4.0.3.csh
endif

# ParMGridGen library if available
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if ( $?PARMGRIDGEN_SYSTEM == 0 && $?WM_THIRD_PARTY_USE_PARMGRIDGEN_10 != 0 && -e $WM_THIRD_PARTY_DIR/packages/ParMGridGen-1.0/platforms/$WM_OPTIONS ) then
    _foamSource $WM_THIRD_PARTY_DIR/packages/ParMGridGen-1.0/platforms/$WM_OPTIONS/etc/ParMGridGen-1.0.csh
endif

# Load Libccmio library
# ~~~~~~~~~~~~~~~~~~~~~
if ( $?LIBCCMIO_SYSTEM == 0 && $?WM_THIRD_PARTY_USE_LIBCCMIO_261 != 0 && -e $WM_THIRD_PARTY_DIR/packages/libccmio-2.6.1/platforms/$WM_OPTIONS ) then
    _foamSource $WM_THIRD_PARTY_DIR/packages/libccmio-2.6.1/platforms/$WM_OPTIONS/etc/libccmio-2.6.1.csh
endif


# Scotch library if available
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~
if ( $?SCOTCH_SYSTEM == 0 && $?WM_THIRD_PARTY_USE_SCOTCH_604 != 0 && -e $WM_THIRD_PARTY_DIR/packages/scotch-6.0.4/platforms/$WM_OPTIONS ) then
    _foamSource $WM_THIRD_PARTY_DIR/packages/scotch-6.0.4/platforms/$WM_OPTIONS/etc/scotch-6.0.4.csh
endif

if ( $?SCOTCH_SYSTEM == 0 && $?WM_THIRD_PARTY_USE_SCOTCH_600 != 0 && -e $WM_THIRD_PARTY_DIR/packages/scotch-6.0.0/platforms/$WM_OPTIONS ) then
    _foamSource $WM_THIRD_PARTY_DIR/packages/scotch-6.0.0/platforms/$WM_OPTIONS/etc/scotch-6.0.0.csh
endif

# Switch on the hoard memory allocator if available
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#if ( -f $FOAM_LIBBIN/libhoard.so ) then
#    setenv LD_PRELOAD $FOAM_LIBBIN/libhoard.so:${LD_PRELOAD}
#endif

# Third party packages

# cmake
# ~~~~~
if ( $?CMAKE_SYSTEM == 0 && $?WM_THIRD_PARTY_USE_CMAKE_322 != 0 && -e "$WM_THIRD_PARTY_DIR"/packages/cmake-3.2.2/platforms/$WM_OPTIONS ) then
    _foamSource $WM_THIRD_PARTY_DIR/packages/cmake-3.2.2/platforms/$WM_OPTIONS/etc/cmake-3.2.2.csh
endif

if ( $?CMAKE_SYSTEM == 0 && $?WM_THIRD_PARTY_USE_CMAKE_2812 != 0 && -e "$WM_THIRD_PARTY_DIR"/packages/cmake-2.8.12/platforms/$WM_OPTIONS ) then
    _foamSource $WM_THIRD_PARTY_DIR/packages/cmake-2.8.12/platforms/$WM_OPTIONS/etc/cmake-2.8.12.csh
endif

# m4
# ~~~~~
if ( $?M4_SYSTEM == 0 && $?WM_THIRD_PARTY_USE_M4_1416 != 0 && -e "$WM_THIRD_PARTY_DIR"/packages/m4-1.4.16/platforms/$WM_OPTIONS ) then
    _foamSource $WM_THIRD_PARTY_DIR/packages/m4-1.4.16/platforms/$WM_OPTIONS/etc/m4-1.4.16.csh
endif

# bison
# ~~~~~
if ( $?BISON_SYSTEM == 0 && $?WM_THIRD_PARTY_USE_BISON_27 != 0 && -e "$WM_THIRD_PARTY_DIR"/packages/bison-2.7/platforms/$WM_OPTIONS ) then
    _foamSource $WM_THIRD_PARTY_DIR/packages/bison-2.7/platforms/$WM_OPTIONS/etc/bison-2.7.csh
endif

# flex
# ~~~~~
if ( $?FLEX_SYSTEM == 0 && $?WM_THIRD_PARTY_USE_FLEX_2535 != 0 && -e "$WM_THIRD_PARTY_DIR"/packages/flex-2.5.35/platforms/$WM_OPTIONS ) then
    _foamSource $WM_THIRD_PARTY_DIR/packages/flex-2.5.35/platforms/$WM_OPTIONS/etc/flex-2.5.35.csh
endif

# zoltan
# ~~~~~
if ( $?ZOLTAN_SYSTEM == 0 && $?WM_THIRD_PARTY_USE_ZOLTAN_35 != 0 && -e "$WM_THIRD_PARTY_DIR"/packages/zoltan-3.5/platforms/$WM_OPTIONS ) then
    _foamSource $WM_THIRD_PARTY_DIR/packages/zoltan-3.5/platforms/$WM_OPTIONS/etc/zoltan-3.5.csh
endif

# Python
# ~~~~~
if ( $?PYTHON_SYSTEM == 0 && $?WM_THIRD_PARTY_USE_PYTHON_27 != 0 && -e "$WM_THIRD_PARTY_DIR"/packages/Python-2.7/platforms/$WM_OPTIONS ) then
    _foamSource $WM_THIRD_PARTY_DIR/packages/Python-2.7/platforms/$WM_OPTIONS/etc/Python-2.7.csh
endif

# PyFoam
# ~~~~~~
if ( $?PYFOAM_SYSTEM == 0 && $?WM_THIRD_PARTY_USE_PYFOAM_063 != 0 && -e "$WM_THIRD_PARTY_DIR"/packages/PyFoam-0.6.3/platforms/noarch ) then
    _foamSource $WM_THIRD_PARTY_DIR/packages/PyFoam-0.6.3/platforms/noarch/etc/PyFoam-0.6.3.csh
endif
if ( $?PYFOAM_SYSTEM == 0 && $?WM_THIRD_PARTY_USE_PYFOAM_064 != 0 && -e "$WM_THIRD_PARTY_DIR"/packages/PyFoam-0.6.4/platforms/noarch ) then
    _foamSource $WM_THIRD_PARTY_DIR/packages/PyFoam-0.6.4/platforms/noarch/etc/PyFoam-0.6.4.csh
endif

# hwloc
# ~~~~~
if ( $?HWLOC_SYSTEM == 0 && $?WM_THIRD_PARTY_USE_HWLOC_1101 != 0 && -e "$WM_THIRD_PARTY_DIR"/packages/hwloc-1.10.1/platforms/$WM_OPTIONS ) then
    _foamSource $WM_THIRD_PARTY_DIR/packages/hwloc-1.10.1/platforms/$WM_OPTIONS/etc/hwloc-1.10.1.csh
endif
if ( $?HWLOC_SYSTEM == 0 && $?WM_THIRD_PARTY_USE_HWLOC_172 != 0 && -e "$WM_THIRD_PARTY_DIR"/packages/hwloc-1.7.2/platforms/$WM_OPTIONS ) then
    _foamSource $WM_THIRD_PARTY_DIR/packages/hwloc-1.7.2/platforms/$WM_OPTIONS/etc/hwloc-1.7.2.csh
endif

# QT
# ~~~~~
if ( $?QT_SYSTEM == 0 && $?WM_THIRD_PARTY_USE_QT_486 != 0 && -e "$WM_THIRD_PARTY_DIR"/packages/qt-everywhere-opensource-src-4.8.6/platforms/$WM_OPTIONS )then
    _foamSource $WM_THIRD_PARTY_DIR/packages/qt-everywhere-opensource-src-4.8.6/platforms/$WM_OPTIONS/etc/qt-everywhere-opensource-src-4.8.6.csh
endif

# PARAVIEW
# ~~~~~
if ( $?PARAVIEW_SYSTEM == 0 && $?WM_THIRD_PARTY_USE_PARAVIEW_440 != 0 && -e "$WM_THIRD_PARTY_DIR"/packages/ParaView-4.4.0/platforms/$WM_OPTIONS ) then
    _foamSource $WM_THIRD_PARTY_DIR/packages/ParaView-4.4.0/platforms/$WM_OPTIONS/etc/ParaView-4.4.0.csh
endif
if ( $?PARAVIEW_SYSTEM == 0 && $?WM_THIRD_PARTY_USE_PARAVIEW_431 != 0 && -e "$WM_THIRD_PARTY_DIR"/packages/ParaView-4.3.1/platforms/$WM_OPTIONS ) then
    _foamSource $WM_THIRD_PARTY_DIR/packages/ParaView-4.3.1/platforms/$WM_OPTIONS/etc/ParaView-4.3.1.csh
endif
if ( $?PARAVIEW_SYSTEM == 0 && $?WM_THIRD_PARTY_USE_PARAVIEW_410 != 0 && -e "$WM_THIRD_PARTY_DIR"/packages/ParaView-4.1.0/platforms/$WM_OPTIONS ) then
    _foamSource $WM_THIRD_PARTY_DIR/packages/ParaView-4.1.0/platforms/$WM_OPTIONS/etc/ParaView-4.1.0.csh
endif
if ( $?PARAVIEW_SYSTEM == 0 && $?WM_THIRD_PARTY_USE_PARAVIEW_401 != 0 && -e "$WM_THIRD_PARTY_DIR"/packages/ParaView-4.0.1/platforms/$WM_OPTIONS ) then
    _foamSource $WM_THIRD_PARTY_DIR/packages/ParaView-4.0.1/platforms/$WM_OPTIONS/etc/ParaView-4.0.1.csh
endif
if ( $?PARAVIEW_SYSTEM == 0 && $?WM_THIRD_PARTY_USE_PARAVIEW_381 != 0 && -e "$WM_THIRD_PARTY_DIR"/packages/ParaView-3.8.1/platforms/$WM_OPTIONS ) then
    _foamSource $WM_THIRD_PARTY_DIR/packages/ParaView-3.8.1/platforms/$WM_OPTIONS/etc/ParaView-3.8.1.csh
endif
if ( $?PARAVIEW_SYSTEM == 0 && $?WM_THIRD_PARTY_USE_PARAVIEW_3141 != 0 && -e "$WM_THIRD_PARTY_DIR"/packages/ParaView-3.14.1/platforms/$WM_OPTIONS ) then
   _foamSource $WM_THIRD_PARTY_DIR/packages/ParaView-3.14.1/platforms/$WM_OPTIONS/etc/ParaView-3.14.1.csh
endif
if ( $?PARAVIEW_SYSTEM == 0 && $?WM_THIRD_PARTY_USE_PARAVIEW_381 != 0 && -e "$WM_THIRD_PARTY_DIR"/packages/ParaView-3.8.1/platforms/$WM_OPTIONS ) then
    _foamSource $WM_THIRD_PARTY_DIR/packages/ParaView-3.8.1/platforms/$WM_OPTIONS/etc/ParaView-3.8.1.csh
endif


# llvm
# ~~~~~
if ( $?LLVM_SYSTEM == 0 && $?WM_THIRD_PARTY_USE_LLVM_360 != 0 && -e "$WM_THIRD_PARTY_DIR"/packages/llvm-3.6.0/platforms/$WM_OPTIONS ) then
    _foamSource $WM_THIRD_PARTY_DIR/packages/llvm-3.6.0/platforms/$WM_OPTIONS/etc/llvm-3.6.0.csh
endif

# mesa
# ~~~~~
if ( $?MESA_SYSTEM == 0 && $?WM_THIRD_PARTY_USE_MESA_GIT != 0 && -e "$WM_THIRD_PARTY_DIR"/packages/mesa-git/platforms/$WM_OPTIONS ) then
    _foamSource $WM_THIRD_PARTY_DIR/packages/mesa-git/platforms/$WM_OPTIONS/etc/mesa-git.csh
endif


if ( $WM_ARCH == "darwinIntel" || $WM_ARCH == "darwinIntel64" ) then
    setenv DYLD_LIBRARY_PATH ${LD_LIBRARY_PATH}
endif

# cleanup environment:
# ~~~~~~~~~~~~~~~~~~~~
unalias _foamAddPath
unalias _foamAddLib
unset minBufferSize

# -----------------------------------------------------------------------------
