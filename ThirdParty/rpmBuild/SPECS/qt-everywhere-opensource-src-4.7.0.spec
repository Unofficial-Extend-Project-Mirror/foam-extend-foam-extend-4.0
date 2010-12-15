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
#     RPM spec file for qt-everywhere-opensource-src-4.7.0
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

%define name		qt-everywhere-opensource-src
%define release		%{_WM_OPTIONS}
%define version 	4.7.0

%define buildroot       %{_topdir}/BUILD/%{name}-%{version}-root

%define _unpackaged_files_terminate_build  0
%define _missing_doc_files_terminate_build 0

BuildRoot:	        %{buildroot}
Summary: 		qt-everywhere-opensource-src
License: 		Unkown
Name: 			%{name}
Version: 		%{version}
Release: 		%{release}
URL:                    http://get.qt.nokia.com/qt/source/qt-everywhere-opensource-src
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

    ./configure                           \
        -opensource --confirm-license=yes \
        --prefix=%{_installPrefix}

    [ -z "$WM_NCOMPPROCS" ] && WM_NCOMPPROCS=1
    make -j $WM_NCOMPPROCS

%install
#    make install INSTALL_ROOT=$RPM_BUILD_ROOT%{_installPrefix}
    make install INSTALL_ROOT=$RPM_BUILD_ROOT

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

export QT_DIR=\$WM_THIRD_PARTY_DIR/packages/%{name}-%{version}/platforms/\$WM_OPTIONS
export QT_BIN_DIR=\$QT_DIR/bin

# Enable access to the runtime package applications
[ -d \$QT_BIN_DIR ] && _foamAddPath \$QT_BIN_DIR
DOT_SH_EOF

    #
    # Generate package specific .csh file for OpenFOAM
    #
cat << DOT_CSH_EOF > $RPM_BUILD_ROOT/%{_installPrefix}/etc/%{name}-%{version}.csh
# Load %{name}-%{version} libraries and binaries if available
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setenv QT_DIR \$WM_THIRD_PARTY_DIR/packages/%{name}-%{version}/platforms/\$WM_OPTIONS
setenv QT_BIN_DIR \$QT_DIR/bin

if ( -e \$QT_BIN_DIR ) then
    _foamAddPath \$QT_BIN_DIR
endif
DOT_CSH_EOF

%clean
rm -rf %{buildroot}

%files
%defattr(-,root,root)
%{_installPrefix}/bin
%{_installPrefix}/etc
%{_installPrefix}/lib
%{_installPrefix}/imports
%{_installPrefix}/include
%{_installPrefix}/mkspecs
%{_installPrefix}/plugins
%{_installPrefix}/phrasebooks
%{_installPrefix}/q3porting.xml
%{_installPrefix}/translations


