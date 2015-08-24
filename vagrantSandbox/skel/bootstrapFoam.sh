#! /usr/bin/env bash

# Just to be sure
export WM_SCHEDULER=ccache
export CCACHE_DIR=/vagrant/ccache4vm

BOOTSTRAPLOG=/home/vagrant/bootstrapFoam.log

cd foam/foam-extend-3.2
source etc/bashrc

( cd wmake/src && make )
cd $WM_THIRD_PARTY_DIR

./AllMake.stage0  2>&1 | tee $BOOTSTRAPLOG
./AllMake.stage1  2>&1 | tee --append $BOOTSTRAPLOG
./AllMake.stage2  2>&1 | tee --append $BOOTSTRAPLOG
source $WM_PROJECT_DIR/etc/bashrc
if [ ! -e $MPI_ARCH_PATH/lib ]
then
    # OpenSUSE needs this
    ln -s $MPI_ARCH_PATH/lib64 $MPI_ARCH_PATH/lib
fi
./AllMake.stage3  2>&1 | tee --append $BOOTSTRAPLOG

cd $WM_PROJECT_DIR
# pick up installed packages
source etc/bashrc

./Allwmake  2>&1 | tee --append $BOOTSTRAPLOG

# compile swak4Foam
cd $WM_THIRD_PARTY_DIR
./AllMake.stage5  2>&1 | tee --append $BOOTSTRAPLOG

# compile the Bazaar
cd $WM_PROJECT_DIR/extend-bazaar
./Allwmake  2>&1 | tee --append $BOOTSTRAPLOG
