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

    # Download prequisite libraries
    ./contrib/download_prerequisites

    # Fixing an old bug related to gmp and m4 and gcc
    # See https://gmplib.org/list-archives/gmp-bugs/2008-August/001114.html
    # This bug will occur when compiling gmp as an internal gcc package.
    # This bug will occur on Centos 6.2. and probably on other OSes as well
    cp gmp/configure.in gmp/configure.in.orig
    cp gmp/configure gmp/configure.orig
    sed "s/m4-not-needed/m4/g" gmp/configure.in > gmp/configure.in.patched 
    sed "s/m4-not-needed/m4/g" gmp/configure > gmp/configure.patched 
    mv gmp/configure.in.patched gmp/configure.in
    mv gmp/configure.patched gmp/configure

    # Create working directory
    mkdir ./objBuildDir
    cd ./objBuildDir

%ifos darwin
    # Couple of fixes for OSX Yosemite and Xcode 6.3.1 
    # See here : http://lists.gnu.org/archive/html/libtool-patches/2014-09/msg00001.html
    # and here : https://gcc.gnu.org/bugzilla/show_bug.cgi?id=63610
    # for the details.
    # Hopefully, this will all be resolved soon before MACOSX_DEPLOYMENT_TARGET=10.11...
    if [ "$MACOSX_DEPLOYMENT_TARGET" = "10.10" ]
    then
	export MACOSX_DEPLOYMENT_TARGET="10.9"
    fi

    # Another necessary fix:
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
