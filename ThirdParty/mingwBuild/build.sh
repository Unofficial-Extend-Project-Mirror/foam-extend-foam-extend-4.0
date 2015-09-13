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
#     build.sh
#
# Description
#     Downloads, extracts, builds and installs thirdy-party dependencies.
#
# Author:
#     Cesare Guardino, Alstom Power Ltd., (2015)
#
#------------------------------------------------------------------------------

# {{{ DEFINE UTILITY FUNCTIONS
download() {
    file=$1
    url=$2

    if [ -f $BUILD_HOME/downloads/$file ] ; then
        echo "Using already existing file $BUILD_HOME/downloads/$file"
    else
        wget --no-check-certificate $url -O $BUILD_HOME/downloads/$file
    fi
}

extract() {
    file=$1
    program=$2

    cp -p $BUILD_HOME/downloads/$file .
    package=`basename $file`
    if [ "$program" = "7zip" ] ; then
        "$ZIP_EXE" x $package
    else
        $program -cd $package | tar xvf -
    fi
    rm $package
}

unzip_dir() {
    dir=$1

    mkdir $dir
    cd $dir
    extract $dir.zip 7zip
    cd ..
}

patch() {
    dir=$1

    cp -rp $BUILD_HOME/$ARCH/patches/$dir .
}

mkchk() {
    dir=$1

    if [ ! -d $dir ] ; then
        mkdir $dir
    fi
}

mkdel() {
    dir=$1

    rm -rf $dir > /dev/null 2>&1
    mkdir $dir
}
# }}}

# {{{ DEFINE PROCESS FUNCTIONS
start() {
    echo "======================== FOAM-EXTEND THIRD-PARTY DEPENDENCIES WINDOWS BUILD SCRIPT ========================"
}

initialise() {
    echo ""

    if [ ! "$MINGW_HOME" ] ; then
        echo "*** ERROR: MINGW_HOME environment variable not specified."
        exit 1
    else
        echo "Using MINGW_HOME=$MINGW_HOME"
    fi

    BUILD_HOME=`pwd`
    ZIP_EXE="7z.exe"
    ARCH="x64"

    BUILD_DIR=$BUILD_HOME/$ARCH/build
    INSTALL_DIR=$BUILD_HOME/$ARCH/install
    OUT_DIR=$BUILD_HOME/$ARCH/output

    mkchk $BUILD_HOME/downloads

    echo ""
    echo "All stdout/stderr output is redirected to the directory $OUT_DIR"
    echo "All builds occur in the directory $BUILD_DIR"
    echo "The script will install the completed builds in the directory $INSTALL_DIR"
}

cleanup() {
    echo ""
    echo "Removing previous builds ..."

    mkdel $BUILD_DIR
    mkdel $INSTALL_DIR
    mkdel $OUT_DIR
}

build_library() {
    PACKAGE=$1

    echo "- Building $PACKAGE ..."
    LOG_FILE=$OUT_DIR/$PACKAGE.log
    cd $BUILD_DIR

    case $PACKAGE in

        dlfcn-win32-master)
            download $PACKAGE.zip https://github.com/dlfcn-win32/dlfcn-win32/archive/master.zip > $LOG_FILE 2>&1
            extract $PACKAGE.zip 7zip >> $LOG_FILE 2>&1
            cd $PACKAGE
            ./configure --prefix=$INSTALL_DIR/system >> $LOG_FILE 2>&1
            make >> $LOG_FILE 2>&1
            mkdir $INSTALL_DIR/system
            make install >> $LOG_FILE 2>&1
            ;;

        system)
            cd $INSTALL_DIR
            patch system
            ;;

        pthreads-w32-2-9-1-release)
            download $PACKAGE.zip ftp://sourceware.org/pub/pthreads-win32/pthreads-w32-2-9-1-release.zip > $LOG_FILE 2>&1
            unzip_dir $PACKAGE >> $LOG_FILE 2>&1
            patch $PACKAGE
            ;;

        metis-5.1.0)
            download $PACKAGE.tar.gz http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/$PACKAGE.tar.gz > $LOG_FILE 2>&1
            extract "$PACKAGE.tar.gz" gzip >> $LOG_FILE 2>&1
            patch $PACKAGE
            cd $PACKAGE
            mkdir build/windows
            cd build/windows
            cmake -G "MSYS Makefiles" -DCMAKE_CONFIGURATION-TYPES="Release" -DGKLIB_PATH="../../GKlib" ../.. >> $LOG_FILE 2>&1
            make >> $LOG_FILE 2>&1
            mkdir $INSTALL_DIR/$PACKAGE
            mkdir $INSTALL_DIR/$PACKAGE/bin
            mkdir $INSTALL_DIR/$PACKAGE/include
            mkdir $INSTALL_DIR/$PACKAGE/lib
            cp -p programs/*.exe $INSTALL_DIR/$PACKAGE/bin
            cp -p ../../include/metis.h $INSTALL_DIR/$PACKAGE/include
            cp -p libmetis/libmetis.a $INSTALL_DIR/$PACKAGE/lib
            ;;

        parmetis-4.0.3)
            download $PACKAGE.tar.gz http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/$PACKAGE.tar.gz > $LOG_FILE 2>&1
            extract "$PACKAGE.tar.gz" gzip >> $LOG_FILE 2>&1
            patch $PACKAGE
            cd $PACKAGE
            mkdir build/windows
            cd build/windows
            cmake -G "MSYS Makefiles" -DCMAKE_CONFIGURATION-TYPES="Release" -DGKLIB_PATH="../../metis/GKlib" ../.. >> $LOG_FILE 2>&1
            $BUILD_HOME/parmetis_includes_hack.pl
            make >> $LOG_FILE 2>&1
            mkdir $INSTALL_DIR/$PACKAGE
            mkdir $INSTALL_DIR/$PACKAGE/bin
            mkdir $INSTALL_DIR/$PACKAGE/include
            mkdir $INSTALL_DIR/$PACKAGE/lib
            cp -p programs/*.exe $INSTALL_DIR/$PACKAGE/bin
            cp -p ../../metis/include/metis.h $INSTALL_DIR/$PACKAGE/include
            cp -p ../../include/parmetis.h $INSTALL_DIR/$PACKAGE/include
            cp -p libmetis/libmetis.a $INSTALL_DIR/$PACKAGE/lib
            cp -p libparmetis/libparmetis.a $INSTALL_DIR/$PACKAGE/lib
            ;;

        ParMGridGen-1.0)
            export EXTRA_SYSTEM_HOME=$INSTALL_DIR/system
            download $PACKAGE.tar.gz http://www.mgnet.org/mgnet/Codes/parmgridgen/$PACKAGE.tar.gz > $LOG_FILE 2>&1
            extract "$PACKAGE.tar.gz" gzip >> $LOG_FILE 2>&1
            patch $PACKAGE
            cd $PACKAGE
            make serial   >> $LOG_FILE 2>&1
            make parallel >> $LOG_FILE 2>&1
            mkdir $INSTALL_DIR/$PACKAGE
            mkdir $INSTALL_DIR/$PACKAGE/bin
            mkdir $INSTALL_DIR/$PACKAGE/include
            mkdir $INSTALL_DIR/$PACKAGE/lib
            cp -p *.exe $INSTALL_DIR/$PACKAGE/bin
            cp -p libmgrid.a $INSTALL_DIR/$PACKAGE/lib
            cp -p libparmgrid.a $INSTALL_DIR/$PACKAGE/lib
            cp -p MGridGen/IMlib/libIMlib.a $INSTALL_DIR/$PACKAGE/lib
            cp -p ParMGridGen/IMParMetis-2.0/libIMparmetis.a $INSTALL_DIR/$PACKAGE/lib
            cp -p MGridGen/IMlib/*.h $INSTALL_DIR/$PACKAGE/include
            cp -p MGridGen/Lib/*.h $INSTALL_DIR/$PACKAGE/include
            export EXTRA_SYSTEM_HOME=
            ;;

        scotch_6.0.4)
            export PTHREADS_HOME=$BUILD_DIR/pthreads-w32-2-9-1-release
            download $PACKAGE.tar.gz https://gforge.inria.fr/frs/download.php/34618 > $LOG_FILE 2>&1
            extract "$PACKAGE.tar.gz" gzip >> $LOG_FILE 2>&1
            patch $PACKAGE
            cd $PACKAGE/src
            make scotch   >> $LOG_FILE 2>&1
            make ptscotch >> $LOG_FILE 2>&1
            mkdir $INSTALL_DIR/$PACKAGE
            make install prefix=$INSTALL_DIR/$PACKAGE >> $PACKAGE.log 2>&1
            export PTHREADS_HOME=
            ;;

        mesquite-2.1.2)
            export CPPFLAGS=-fpermissive
            download $PACKAGE.tar.gz http://downloads.sourceforge.net/project/openfoam-extend/foam-extend-3.0/ThirdParty/$PACKAGE.tar.gz > $LOG_FILE 2>&1
            extract "$PACKAGE.tar.gz" gzip >> $LOG_FILE 2>&1
            cd $PACKAGE
            cp -p $MINGW_HOME/bin/libstdc++-6.dll utils
            ./configure --prefix=$INSTALL_DIR >> $LOG_FILE 2>&1
            make >> $LOG_FILE 2>&1
            mkdir $INSTALL_DIR/$PACKAGE
            make install prefix=$INSTALL_DIR/$PACKAGE >> $LOG_FILE 2>&1
            export CPPFLAGS=
            ;;

        *)
            echo "*** ERROR: Unknown package '$PACKAGE'"
            exit 1
            ;;
    esac
}

build_libraries() {
    echo ""
    echo "Building libraries ..."
    build_library dlfcn-win32-master
    build_library system
    build_library pthreads-w32-2-9-1-release
    build_library metis-5.1.0
    build_library parmetis-4.0.3
    build_library ParMGridGen-1.0
    build_library scotch_6.0.4
    build_library mesquite-2.1.2
}

create_dirs() {
    echo ""
    echo "Checking for build directories and creating them if required ..."

    mkchk $BUILD_DIR
    mkchk $INSTALL_DIR
    mkchk $OUT_DIR
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
mkchk $ARCH
cleanup
build_libraries
finish
# }}}
