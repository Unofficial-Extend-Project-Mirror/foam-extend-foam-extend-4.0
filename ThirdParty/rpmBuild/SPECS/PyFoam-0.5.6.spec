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
#     RPM spec file for PyFoam-0.5.6
#
# Description
#     RPM spec file for PyFoam version 0.5.6
#
# Author:
#     Martin Beaudoin, Hydro-Quebec, (2012)
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

%define name		PyFoam
%define release		1
%define version 	0.5.6

%define buildroot       %{_topdir}/BUILD/%{name}-%{version}-root

BuildRoot:	        %{buildroot}
Summary: 		PyFoam
License: 		Unkown
Name: 			%{name}
Version: 		%{version}
Release: 		%{release}
BuildArch:              noarch
URL:                    http://openfoamwiki.net/images/b/b8/PyFoam-0.5.6.tar.gz
Source: 		%url/%{name}-%{version}.tar.gz
Prefix: 		%{_prefix}
Group: 			Development/Tools

%define _installPrefix  %{_prefix}/packages/%{name}-%{version}/platforms/noarch

%description
%{summary}

%prep
%setup -q

%build
    # Nothing to do
    # true

%install
    python setup.py install --prefix=$RPM_BUILD_ROOT/%{_installPrefix}

    %define pythonVersion $(python -V 2>&1 | awk -F ' ' {'print $2'} |  awk -F \. {'print $1 "." $2'})

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

export PYFOAM_DIR=\$WM_THIRD_PARTY_DIR/packages/%{name}-%{version}/platforms/noarch
export PYTHONPATH=\$PYFOAM_DIR/lib/python%{pythonVersion}/site-packages:\$PYTHONPATH

# Enable access to the package applications if present
[ -d \$PYFOAM_DIR/bin ] && _foamAddPath \$PYFOAM_DIR/bin

DOT_SH_EOF

    #
    # Generate package specific .csh file for OpenFOAM
    #
cat << DOT_CSH_EOF > $RPM_BUILD_ROOT/%{_installPrefix}/etc/%{name}-%{version}.csh
# Load %{name}-%{version} libraries and binaries if available
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setenv PYFOAM_DIR \$WM_THIRD_PARTY_DIR/packages/%{name}-%{version}/platforms/noarch

setenv PYTHONPATH \$PYFOAM_DIR/lib/python%{pythonVersion}/site-packages:\$PYTHONPATH

if ( -e \$PYFOAM_DIR/bin ) then
    _foamAddPath \$PYFOAM_DIR/bin
endif

if ( -e \$PYFOAM_SITE_DIR/bin ) then
    _foamAddPath \$PYFOAM_SITE_DIR/bin
endif

DOT_CSH_EOF

cat << DOT_HARDCODED_SH_EOF > $RPM_BUILD_ROOT/%{_installPrefix}/etc/%{name}-%{version}_hardcoded.sh
# In this version of the configuration script, the paths are hardcoded,
# which makes it easier to load PyFoam without the OpenFOAM environment
# variables
#
# Load %{name}-%{version} libraries and binaries if available
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

export PYFOAM_DIR=$WM_THIRD_PARTY_DIR/packages/%{name}-%{version}/platforms/noarch
export PYTHONPATH=\$PYFOAM_DIR/lib/python%{pythonVersion}/site-packages:\$PYTHONPATH

# Enable access to the package applications if present
[ -d \$PYFOAM_DIR/bin ] && export PATH=\$PYFOAM_DIR/bin:\$PATH

DOT_HARDCODED_SH_EOF

    #
    # Generate package specific .csh file for OpenFOAM
    #
cat << DOT_HARDCODED_CSH_EOF > $RPM_BUILD_ROOT/%{_installPrefix}/etc/%{name}-%{version}_hardcoded.csh
# In this version of the configuration script, the paths are hardcoded,
# which makes it easier to load PyFoam without the OpenFOAM environment
# variables
#
# Load %{name}-%{version} libraries and binaries if available
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setenv PYFOAM_DIR $WM_THIRD_PARTY_DIR/packages/%{name}-%{version}/platforms/noarch

setenv PYTHONPATH \$PYFOAM_DIR/lib/python%{pythonVersion}/site-packages:\$PYTHONPATH

if ( -e \$PYFOAM_DIR/bin ) then
    setenv PATH \$PYFOAM_DIR/bin:\$PATH
endif

DOT_HARDCODED_CSH_EOF

    #finally, generate a .tgz file for systems where using rpm for installing packages
    # as a non-root user might be a problem.
    (mkdir -p  %{_topdir}/TGZS/%{_target_cpu}; cd $RPM_BUILD_ROOT/%{_prefix}; tar -zcvf %{_topdir}/TGZS/%{_target_cpu}/%{name}-%{version}.tgz  packages/%{name}-%{version})
 

%clean
rm -rf %{buildroot}

%files
%defattr(-,root,root)
%{_installPrefix}

