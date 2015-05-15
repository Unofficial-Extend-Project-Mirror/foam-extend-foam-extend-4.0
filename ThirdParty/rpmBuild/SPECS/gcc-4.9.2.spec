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
#     RPM spec file for gcc-4.9.2
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
#   By using this prefix for the Prefix:  parameter in thi file, you will make this
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

%define name		gcc
%define release		%{_WM_OPTIONS}
%define version 	4.9.2

%define buildroot       %{_topdir}/BUILD/%{name}-%{version}-root

BuildRoot:	        %{buildroot}
Summary: 		gcc
License: 		Unkown
Name: 			%{name}
Version: 		%{version}
Release: 		%{release}
URL:                    ftp://ftp.gnu.org/gnu/gcc/gcc-4.9.2
Source: 		%url/%{name}-%{version}.tar.gz
Prefix: 		%{_prefix}
Group: 			Development/Tools


%define _installPrefix  %{_prefix}/packages/%{name}-%{version}/platforms/%{_WM_OPTIONS}

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

    GMP_VERSION=gmp-6.0.0
    MPFR_VERSION=mpfr-3.1.2
    MPC_VERSION=mpc-1.0.3
    ISL_VERSION=isl-0.12.2
    CLOOG_VERSION=cloog-0.18.1

    # Download dependency packages and add them as builtin components for gcc
    #GMP
    wget ftp://ftp.gnu.org/gnu/gmp/${GMP_VERSION}a.tar.xz
    tar -xf ${GMP_VERSION}a.tar.xz
    mv ${GMP_VERSION} gmp

    #MPFR
    wget ftp://ftp.gnu.org/gnu/mpfr/$MPFR_VERSION.tar.gz
    tar -xf $MPFR_VERSION.tar.gz
    mv ${MPFR_VERSION} mpfr

    #MPC
    wget ftp://ftp.gnu.org/gnu/mpc/$MPC_VERSION.tar.gz
    tar -xf $MPC_VERSION.tar.gz
    mv ${MPC_VERSION} mpc
    
    #ISL
    wget ftp://gcc.gnu.org/pub/gcc/infrastructure/$ISL_VERSION.tar.bz2
    tar -xf $ISL_VERSION.tar.bz2
    mv ${ISL_VERSION} isl

    #CLOOG
    wget ftp://gcc.gnu.org/pub/gcc/infrastructure/$CLOOG_VERSION.tar.gz
    tar -xf $CLOOG_VERSION.tar.gz
    mv ${CLOOG_VERSION} cloog

    mkdir ./objBuildDir
    cd ./objBuildDir

%ifos darwin
    # Use 'bootstrap-debug' build configuration to force stripping of object
    # files prior to comparison during bootstrap (broken by Xcode 6.3).
    # Fix taken from HomeBrew and MacPorts projects
    ../configure     \
        --prefix=%{_installPrefix} \
        --enable-languages=c,c++   \
        --enable-shared            \
        --with-build-config=bootstrap-debug \
	--disable-multilib
%else
    ../configure     \
        --prefix=%{_installPrefix} \
        --enable-languages=c,c++   \
        --enable-shared            \
	--disable-multilib
%endif

    [ -z "$WM_NCOMPPROCS" ] && WM_NCOMPPROCS=1
    make -j $WM_NCOMPPROCS

%install
    cd ./objBuildDir
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
# Load %{name}-%{version} libraries and binaries if available
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

export GCC_DIR=\$WM_THIRD_PARTY_DIR/packages/%{name}-%{version}/platforms/\$WM_OPTIONS

[ -d \$GCC_DIR/lib ] && _foamAddLib \$GCC_DIR/lib

[ -d \$GCC_DIR/lib64 ] && _foamAddLib \$GCC_DIR/lib64

# Enable access to the package applications if present
[ -d \$GCC_DIR/bin ] && _foamAddPath \$GCC_DIR/bin
DOT_SH_EOF

    #
    # Generate package specific .csh file for foam-extend
    #
cat << DOT_CSH_EOF > $RPM_BUILD_ROOT/%{_installPrefix}/etc/%{name}-%{version}.csh
# Load %{name}-%{version} libraries and binaries if available
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setenv GCC_DIR \$WM_THIRD_PARTY_DIR/packages/%{name}-%{version}/platforms/\$WM_OPTIONS

if ( -e \$GCC_DIR/lib ) then
    _foamAddLib \$GCC_DIR/lib
endif

if ( -e \$GCC_DIR/lib64 ) then
    _foamAddLib \$GCC_DIR/lib64
endif

if ( -e \$GCC_DIR/bin ) then
    _foamAddPath \$GCC_DIR/bin
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
