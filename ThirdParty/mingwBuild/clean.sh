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
#     clean.sh
#
# Description
#     Removes previous thirdy-party dependencies build directories
#     (does not remove installed packages directory)
#
# Author:
#     Cesare Guardino, Alstom Power Ltd., (2015)
#
#------------------------------------------------------------------------------

# {{{ DEFINE UTILITY FUNCTIONS
remove_dir() {
    dir=$1

    rm -rf $dir > /dev/null 2>&1
}
# }}}

# {{{ DEFINE PROCESS FUNCTIONS
start() {
    echo "======================== FOAM-EXTEND THIRD-PARTY DEPENDENCIES WINDOWS CLEAN SCRIPT ========================"
}

initialise() {
    echo ""

    BUILD_HOME=`pwd`
    ARCH="x64"

    BUILD_DIR=$BUILD_HOME/$ARCH/build
    INSTALL_DIR=$BUILD_HOME/$ARCH/install
    OUT_DIR=$BUILD_HOME/$ARCH/output
}

cleanup() {
    echo ""
    echo "Removing previous builds ..."

    remove_dir $BUILD_DIR
    remove_dir $INSTALL_DIR
    remove_dir $OUT_DIR
    remove_dir $BUILD_HOME/downloads
}

finish() {
    echo ""
    echo "All done!"
}
# }}}

# {{{ MAIN EXECUTION
cd ${0%/*} || exit 1    # run from this directory
start
initialise
cleanup
finish
# }}}
