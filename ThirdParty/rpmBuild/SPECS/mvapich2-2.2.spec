#------------------------------------------------------------------------------
# =========                 |
# \\      /  F ield         | foam-extend: Open Source CFD
#  \\    /   O peration     |
#   \\  /    A nd           | For copyright notice see file Copyright
#    \\/     M anipulation  |
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
#     RPM spec file for mvapich2-2.2
#
# Description
#     RPM spec file for creating a relocatable RPM
#
# Author:
#     Martin Beaudoin, Hydro-Quebec, (2015)
#
#------------------------------------------------------------------------------

# We grab the value of WM_THIRD_PARTY and WM_OPTIONS from the environment variable
%{expand:%%define _WM_THIRD_PARTY_DIR %(echo $WM_THIRD_PARTY_DIR)}
%{expand:%%define _WM_OPTIONS         %(echo $WM_OPTIONS)}

# Disable the generation of debuginfo packages
%define debug_package %{nil}

# The topdir needs to point to the $WM_THIRD_PARTY/rpmbuild directory
%define _topdir	 	%{_WM_THIRD_PARTY_DIR}/rpmBuild
%define _tmppath	%{_topdir}/tmp

# Will install the package directly $WM_THIRD_PARTY_DIR
#   Some comments about package relocation:
#   By using this prefix for the Prefix:  parameter in this file, you will make this
#   package relocatable.
#
#   This is fine, as long as your software is itself relocatable.
#
#   Simply take note that libraries built with libtool are not relocatable because the
#   prefix we specify will be hard-coded in the library .la files.
#   Ref: http://sourceware.org/autobook/autobook/autobook_80.html
#
#   In that case, if you ever change the value of the $WM_THIRD_PARTY_DIR, you will
#   not be able to reutilize this RPM, even though it is relocatable. You will need to
#   regenerate the RPM.
#
%define _prefix         %{_WM_THIRD_PARTY_DIR}

%define name		mvapich2
%define release		%{_WM_OPTIONS}
%define version 	2.2

%define buildroot       %{_topdir}/BUILD/%{name}-%{version}-root

BuildRoot:	        %{buildroot}
Summary: 		mvapich2
License: 		BSD
Name: 			%{name}
Version: 		%{version}
Release: 		%{release}
URL:            http://mvapich.cse.ohio-state.edu/download/mvapich/mv2/
Source: 		%url/%{name}-%{version}.tar.gz
Prefix: 		%{_prefix}
Group: 			Development/Tools


%define _installPrefix  %{_prefix}/packages/%{name}-%{version}/platforms/%{_WM_OPTIONS}

#--------------------------------------------------------------------------
#
# Here, we define default compiling options for mvapich2
#
# One can override the option on the commande line : --define='MACRO EXPR'
#
%{!?_configureAdditionalArgs: %define _configureAdditionalArgs Undefined}


%description
%{summary}

%prep
%setup -q


%build
    # export WM settings in a form that GNU configure recognizes
    [ -n "$WM_CC" ]         &&  export CC="$WM_CC"
    [ -n "$WM_CXX" ]        &&  export CXX="$WM_CXX"
    [ -n "$WM_CFLAGS" ]     &&  export CFLAGS="$WM_CFLAGS"
    [ -n "$WM_CXXFLAGS" ]   &&  export CXXFLAGS="$WM_CXXFLAGS"
    [ -n "$WM_LDFLAGS" ]    &&  export LDFLAGS="$WM_LDFLAGS"

    set +x
    echo ""
    echo "Compilation options:"
    echo "     _configureAdditionalArgs : %{_configureAdditionalArgs}"
    echo ""
    set -x

    unset mpiWith
    # Enable GridEngine if it appears to be in use
    # If you don't want any integration with SGE, simply unset the SGE
    # environment variable
    if [ -n "$SGE_ROOT" ]
    then
        mpiWith="$mpiWith --with-sge"
    else
        mpiWith="$mpiWith --without-sge"
	mpiWith="$mpiWith --enable-mca-no-build=ras-gridengine,pls-gridengine"
    fi

    # Infiniband support
    # Adjust according to your local paths
    # if [ -d /usr/local/ofed -a -d /usr/local/ofed/lib64 ]
    # then
    #     mpiWith="$mpiWith --with-openib=/usr/local/ofed"
    #     mpiWith="$mpiWith --with-openib-libdir=/usr/local/ofed/lib64"
    # fi

    # NOTE. Option --disable-mcast disables Infiniband when libs are not present

    ./configure \
        --prefix=%{_installPrefix}  \
        --exec_prefix=%{_installPrefix}  \
        --enable-mpirun-prefix-by-default \
        --enable-orterun-prefix-by-default \
        --enable-shared \
        --enable-mpi-cxx \
        --disable-static \
        --disable-mpi-f77 \
        --disable-mpi-f90 \
        --without-slurm \
        --enable-mpi-profile $mpiWith \
        --disable-vt \
        echo %{?_configureAdditionalArgs}`

    [ -z "$WM_NCOMPPROCS" ] && WM_NCOMPPROCS=1
    make -j $WM_NCOMPPROCS

%install
    make install DESTDIR=$RPM_BUILD_ROOT

    # Creation of foam-extend specific .csh and .sh files"

    echo ""
    echo "Generating foam-extend specific .csh and .sh files for the package %{name}-%{version}"
    echo ""
    #
    # Generate package specific .sh file for foam-extend
    #
mkdir -p $RPM_BUILD_ROOT/%{_installPrefix}/etc
cat << DOT_SH_EOF > $RPM_BUILD_ROOT/%{_installPrefix}/etc/%{name}-%{version}.sh
# Load %{name}-%{version} libraries and binaries
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
export MVAPICH2_DIR=\$WM_THIRD_PARTY_DIR/packages/%{name}-%{version}/platforms/\$WM_OPTIONS
export MVAPICH2_BIN_DIR=\$MVAPICH2_DIR/bin
export MVAPICH2_LIB_DIR=\$MVAPICH2_DIR/lib

# Enable access to the package runtime applications and libraries
[ -d \$MVAPICH2_BIN_DIR ] && _foamAddPath \$MVAPICH2_BIN_DIR
[ -d \$MVAPICH2_LIB_DIR ] && _foamAddLib  \$MVAPICH2_LIB_DIR

export MPI_HOME=\$MVAPICH2_DIR
export MPI_ARCH_PATH=\$MPI_HOME
export OPAL_PREFIX=\$MPI_ARCH_PATH

# We initialize the rest of the environment using mpicc --showme:
# This does not work with MVAPICH2
#export MVAPICH2_INCLUDE_DIR="\`mpicc --showme:incdirs\`"
#export MVAPICH2_COMPILE_FLAGS="\`mpicc --showme:compile\`"
#export MVAPICH2_LINK_FLAGS="\`mpicc --showme:link\`"

# Set the foam-extend compilation flags
export PINC=\$MVAPICH2_COMPILE_FLAGS
export PLIBS=\$MVAPICH2_LINK_FLAGS


if [ "\$FOAM_VERBOSE" -a "\$PS1" ]
then
    echo "  Environment variables defined for OpenMPI:"
    echo "    MVAPICH2_BIN_DIR       : \$MVAPICH2_BIN_DIR"
    echo "    MVAPICH2_LIB_DIR       : \$MVAPICH2_LIB_DIR"
    echo "    MVAPICH2_INCLUDE_DIR   : \$MVAPICH2_INCLUDE_DIR"
    echo "    MVAPICH2_COMPILE_FLAGS : \$MVAPICH2_COMPILE_FLAGS"
    echo "    MVAPICH2_LINK_FLAGS    : \$MVAPICH2_LINK_FLAGS"
    echo ""
    echo "    MPI_HOME              : \$MPI_HOME"
    echo "    MPI_ARCH_PATH         : \$MPI_ARCH_PATH"
    echo "    OPAL_PREFIX           : \$OPAL_PREFIX"
    echo "    PINC                  : \$PINC"
    echo "    PLIBS                 : \$PLIBS"
fi
DOT_SH_EOF

    #
    # Generate package specific .csh file for foam-extend
    #
cat << DOT_CSH_EOF > $RPM_BUILD_ROOT/%{_installPrefix}/etc/%{name}-%{version}.csh
# Load %{name}-%{version} libraries and binaries
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setenv MVAPICH2_DIR \$WM_THIRD_PARTY_DIR/packages/%{name}-%{version}/platforms/\$WM_OPTIONS
setenv MVAPICH2_BIN_DIR \$MVAPICH2_DIR/bin
setenv MVAPICH2_LIB_DIR \$MVAPICH2_DIR/lib

# Enable access to the package runtime applications and libraries
if ( -e \$MVAPICH2_BIN_DIR ) then
    _foamAddPath \$MVAPICH2_BIN_DIR
endif

if ( -e \$MVAPICH2_LIB_DIR ) then
    _foamAddLib \$MVAPICH2_LIB_DIR
endif

setenv MPI_HOME \$MVAPICH2_DIR
setenv MPI_ARCH_PATH \$MPI_HOME
setenv OPAL_PREFIX \$MPI_ARCH_PATH

# We initialize the rest of the environment using mpicc --showme:
setenv MVAPICH2_INCLUDE_DIR "\`mpicc --showme:incdirs\`"
setenv MVAPICH2_COMPILE_FLAGS "\`mpicc --showme:compile\`"
setenv MVAPICH2_LINK_FLAGS "\`mpicc --showme:link\`"

# Set the foam-extend compilation flags
setenv PINC "\$MVAPICH2_COMPILE_FLAGS"
setenv PLIBS "\$MVAPICH2_LINK_FLAGS"


if (\$?FOAM_VERBOSE && \$?prompt) then
    echo "  Environment variables defined for OpenMPI:"
    echo "    MVAPICH2_BIN_DIR       : \$MVAPICH2_BIN_DIR"
    echo "    MVAPICH2_LIB_DIR       : \$MVAPICH2_LIB_DIR"
    echo "    MVAPICH2_INCLUDE_DIR   : \$MVAPICH2_INCLUDE_DIR"
    echo "    MVAPICH2_COMPILE_FLAGS : \$MVAPICH2_COMPILE_FLAGS"
    echo "    MVAPICH2_LINK_FLAGS    : \$MVAPICH2_LINK_FLAGS"
    echo ""
    echo "    MPI_HOME              : \$MPI_HOME"
    echo "    MPI_ARCH_PATH         : \$MPI_ARCH_PATH"
    echo "    OPAL_PREFIX           : \$OPAL_PREFIX"
    echo "    PINC                  : \$PINC"
    echo "    PLIBS                 : \$PLIBS"
endif
DOT_CSH_EOF

    #finally, generate a .tgz file for systems where using rpm for installing packages
    # as a non-root user might be a problem.
    (mkdir -p  %{_topdir}/TGZS/%{_target_cpu}; cd $RPM_BUILD_ROOT/%{_prefix}; tar -zcvf %{_topdir}/TGZS/%{_target_cpu}/%{name}-%{version}.tgz  packages/%{name}-%{version})


%clean
rm -rf %{buildroot}

%files
%defattr(-,root,root)
%{_installPrefix}
