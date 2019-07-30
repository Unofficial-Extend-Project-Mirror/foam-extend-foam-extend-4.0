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
#     RPM spec file for metis-5.1.0
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

%define name		metis
%define release		%{_WM_OPTIONS}
%define version 	5.1.0

%define buildroot       %{_topdir}/BUILD/%{name}-%{version}-root

BuildRoot:	    %{buildroot}
Summary: 		metis
License: 		Unkown
Name: 			%{name}
Version: 		%{version}
Release: 		%{release}
URL:            http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis
Source: 		%url/%{name}-%{version}.tar.gz
Prefix: 		%{_prefix}
Group: 			Development/Tools
Patch0:         metis-5.1.0_patch_gcc
Patch1:         metis-5.1.0_patch_64Bit

%define _installPrefix  %{_prefix}/packages/%{name}-%{version}/platforms/%{_WM_OPTIONS}

%description
%{summary}

%prep
%setup -q

%ifos linux
%patch0 -p1
%endif

if [ "$WM_LABEL_SIZE" = "64" ]; then
%patch1 -p1
fi

%build
    [ -n "$WM_CC" ]         &&  export CC="$WM_CC"
    [ -n "$WM_CXX" ]        &&  export CXX="$WM_CXX"
    [ -n "$WM_CFLAGS" ]     &&  export CFLAGS="$WM_CFLAGS"
    [ -n "$WM_CXXFLAGS" ]   &&  export CXXFLAGS="$WM_CXXFLAGS"
    [ -n "$WM_LDFLAGS" ]    &&  export LDFLAGS="$WM_LDFLAGS"

    [ -n "$WM_CXX" ]        &&  export OMPI_CXX="$WM_CXX"

    [ -z "$WM_NCOMPPROCS" ] && WM_NCOMPPROCS=1

%ifos darwin
    make config
    # The parameter -D_POSIX_C_SOURCE=200809 has a side effect on Mac OS X
    make -j $WM_NCOMPPROCS OPTFLAGS="-O3 -fPIC" CC="$CC"
%else
    make config
    make -j $WM_NCOMPPROCS OPTFLAGS="-O3 -fPIC -D_POSIX_C_SOURCE=200809" CC="$CC"
%endif

%install
    # Manual installation
    mkdir -p $RPM_BUILD_ROOT/%{_installPrefix}/include
    mkdir -p $RPM_BUILD_ROOT/%{_installPrefix}/lib
    cp build/*/libmetis/libmetis.* $RPM_BUILD_ROOT/%{_installPrefix}/lib
    cp include/metis.h           $RPM_BUILD_ROOT/%{_installPrefix}/include

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

export METIS_DIR=\$WM_THIRD_PARTY_DIR/packages/%{name}-%{version}/platforms/\$WM_OPTIONS
export METIS_BIN_DIR=\$METIS_DIR/bin
export METIS_LIB_DIR=\$METIS_DIR/lib
export METIS_INCLUDE_DIR=\$METIS_DIR/include

# Enable access to the package runtime applications and libraries
[ -d \$METIS_BIN_DIR ] && _foamAddPath \$METIS_BIN_DIR
[ -d \$METIS_LIB_DIR ] && _foamAddLib  \$METIS_LIB_DIR

DOT_SH_EOF

    #
    # Generate package specific .csh file for foam-extend
    #
cat << DOT_CSH_EOF > $RPM_BUILD_ROOT/%{_installPrefix}/etc/%{name}-%{version}.csh
# Load %{name}-%{version} libraries and binaries if available
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setenv METIS_DIR \$WM_THIRD_PARTY_DIR/packages/%{name}-%{version}/platforms/\$WM_OPTIONS
setenv METIS_BIN_DIR \$METIS_DIR/bin
setenv METIS_LIB_DIR \$METIS_DIR/lib
setenv METIS_INCLUDE_DIR \$METIS_DIR/include

if ( -e \$METIS_BIN_DIR ) then
    _foamAddPath \$METIS_BIN_DIR
endif

if ( -e \$METIS_LIB_DIR ) then
    _foamAddLib \$METIS_LIB_DIR
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




