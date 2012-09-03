#------------------------------------------------------------------------------
# =========                 |
# \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
#  \\    /   O peration     |
#   \\  /    A nd           | Copyright held by original author
#    \\/     M anipulation  |
#------------------------------------------------------------------------------
# License
#     This file is part of OpenFOAM.
#
#     OpenFOAM is free software; you can redistribute it and/or modify it
#     under the terms of the GNU General Public License as published by the
#     Free Software Foundation; either version 2 of the License, or (at your
#     option) any later version.
#
#     OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
#     ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#     FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#     for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with OpenFOAM; if not, write to the Free Software Foundation,
#     Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
#
# Script
#     RPM spec file for libccmio-2.6.1
#
# Description
#     RPM spec file for creating a relocatable RPM
#
# Author:
#     Martin Beaudoin, Hydro-Quebec, (2010)
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

%define name		libccmio
%define release		%{_WM_OPTIONS}
%define version 	2.6.1

%define buildroot       %{_topdir}/BUILD/%{name}-%{version}-root

BuildRoot:	        %{buildroot}
Summary: 		libccmio
License: 		Unkown
Name: 			%{name}
Version: 		%{version}
Release: 		%{release}
URL:                    https://wci.llnl.gov/codes/visit/3rd_party
Source: 		%url/%{name}-%{version}.tar.gz
Prefix: 		%{_prefix}
Group: 			Development/Tools
Patch0:                 libccmio-2.6.1.patch_0

%define _installPrefix  %{_prefix}/packages/%{name}-%{version}/platforms/%{_WM_OPTIONS}

%description
%{summary}

%prep
%setup -q
%patch0 -p1

%build
    [ -n "$WM_CC" ]         &&  export CC="$WM_CC"
    [ -n "$WM_CXX" ]        &&  export CXX="$WM_CXX"
    [ -n "$WM_CFLAGS" ]     &&  export CFLAGS="$WM_CFLAGS"
    [ -n "$WM_CXXFLAGS" ]   &&  export CXXFLAGS="$WM_CXXFLAGS"
    [ -n "$WM_LDFLAGS" ]    &&  export LDFLAGS="$WM_LDFLAGS"
    [ -z "$WM_NCOMPPROCS" ] && WM_NCOMPPROCS=1

%ifos darwin
    # Missing configuration files for Mac OS X 
    [ ! -d config/i386-apple-darwin10 ] && cp -r config/i386-apple-darwin8 config/i386-apple-darwin10
    [ ! -d config/i386-apple-darwin11 ] && cp -r config/i386-apple-darwin8 config/i386-apple-darwin11
%endif
    # Warning:
    #  1: The name of the ADF library will be renamed to libadf_ccmio since this
    #     is an old version of ADF, and the author modified the source code
    #  2: The name of the CGNS library will be renamed to libcgns_ccmio as well
    #     since this is an old version of CGNS.
    #
    #  This way, the libraries libadf_ccmio and libcgns_ccmio will not get in
    #  conflict with any other packages that might depend on a newer version
    #  of libadf or libcgns 
    #
    unset RELEASE
    unset DEBUG
    unset STATIC
    unset SHARED
    if [ -d libadf ];    then ( cd libadf;    RELEASE=1 SHARED=1 make -f Makefile.qmake all; ) fi
    if [ -d libccmio ];  then ( cd libccmio;  RELEASE=1 SHARED=1 make -f Makefile.qmake all; ) fi

    # We don't need libcgns_ccmio. We keep it here as a reference
    #if [ -d libcgns ];    then cd libcgns;   RELEASE=1 SHARED=1 make -f Makefile.qmake all;  fi

%install
    # Manual installation
    mkdir -p $RPM_BUILD_ROOT/%{_installPrefix}/include/libccmio
    mkdir -p $RPM_BUILD_ROOT/%{_installPrefix}/lib
    libsdir=`find ./lib -name release-shared`
    mv ${libsdir}/* $RPM_BUILD_ROOT/%{_installPrefix}/lib
    cp libccmio/*.h $RPM_BUILD_ROOT/%{_installPrefix}/include/libccmio

    # Creation of OpenFOAM specific .csh and .sh files"

    echo ""
    echo "Generating OpenFOAM specific .csh and .sh files for the package %{name}-%{version}"
    echo ""
    #
    # Generate package specific .sh file for OpenFOAM
    #
mkdir -p $RPM_BUILD_ROOT/%{_installPrefix}/etc
cat << DOT_SH_EOF > $RPM_BUILD_ROOT/%{_installPrefix}/etc/%{name}-%{version}.sh
# Load %{name}-%{version} libraries and binaries if available
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

export LIBCCMIO_DIR=\$WM_THIRD_PARTY_DIR/packages/%{name}-%{version}/platforms/\$WM_OPTIONS
export LIBCCMIO_BIN_DIR=\$LIBCCMIO_DIR/bin
export LIBCCMIO_LIB_DIR=\$LIBCCMIO_DIR/lib
export LIBCCMIO_INCLUDE_DIR=\$LIBCCMIO_DIR/include

# Enable access to the package runtime applications and libraries
[ -d \$LIBCCMIO_BIN_DIR ] && _foamAddPath \$LIBCCMIO_BIN_DIR
[ -d \$LIBCCMIO_LIB_DIR ] && _foamAddLib  \$LIBCCMIO_LIB_DIR

DOT_SH_EOF

    #
    # Generate package specific .csh file for OpenFOAM
    #
cat << DOT_CSH_EOF > $RPM_BUILD_ROOT/%{_installPrefix}/etc/%{name}-%{version}.csh
# Load %{name}-%{version} libraries and binaries if available
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setenv LIBCCMIO_DIR \$WM_THIRD_PARTY_DIR/packages/%{name}-%{version}/platforms/\$WM_OPTIONS
setenv LIBCCMIO_BIN_DIR \$LIBCCMIO_DIR/bin
setenv LIBCCMIO_LIB_DIR \$LIBCCMIO_DIR/lib
setenv LIBCCMIO_INCLUDE_DIR \$LIBCCMIO_DIR/include

if ( -e \$LIBCCMIO_BIN_DIR ) then
    _foamAddPath \$LIBCCMIO_BIN_DIR
endif

if ( -e \$LIBCCMIO_LIB_DIR ) then
    _foamAddLib \$LIBCCMIO_LIB_DIR
endif
DOT_CSH_EOF

    #finally, generate a .tgz file for systems where using rpm for installing packages
    # as a non-root user might be a problem.
    (mkdir -p  %{_topdir}/TGZS/%{_target_cpu}; cd $RPM_BUILD_ROOT/%{_prefix}; tar -zcvf %{_topdir}/TGZS/%{_target_cpu}/%{name}-%{version}.tgz  packages/%{name}-%{version})


%clean
rm -rf %{buildroot}

%Files
%defattr(-,root,root)
%{_installPrefix}/*




