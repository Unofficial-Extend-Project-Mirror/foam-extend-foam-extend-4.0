#!/bin/sh
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
#     configure_OpenMPI.sh
#
# Description
#     Generates static OpenMPI library to enable compilation with MINGW toolchain.
#
# Author:
#     Cesare Guardino, Alstom Power Ltd., (2015)
#
#------------------------------------------------------------------------------

if [ -f "$MPI_ROOTDIR/lib/libmpi.a" ] ; then
    echo "$MPI_ROOTDIR/lib/libmpi.a already exists."
    exit 1
fi

INSTALL=1
if [ $1 ] ; then
    if [ $1 = "--no-install" ] ; then
        INSTALL=0
        echo "*** WARNING: Will not install exported libmpi.a and libmpi.def"
    fi
fi

tmp=`echo $RANDOM$RANDOM$RANDOM`
current_dir=`pwd`
mkdir -p /tmp/$tmp
cd /tmp/$tmp

pexports $MPI_ROOTDIR/bin/libmpi.dll > libmpi.def
dlltool --dllname libmpi.dll --def libmpi.def --output-lib libmpi.a

if [ $INSTALL -eq 1 ] ; then
    mv libmpi.a $MPI_ROOTDIR/lib
    mv libmpi.def $MPI_ROOTDIR/lib
    cd $current_dir
    rm -rf /tmp/$tmp
    echo "Installed exported libraries into $MPI_ROOTDIR/lib"
else
    cd $current_dir
    echo "Exported libraries left in directory /tmp/$tmp"
fi
