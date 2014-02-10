#! /bin/bash

# Just to be sure
export WM_SCHEDULER=ccache
export CCACHE_DIR=/vagrant/ccache4vm

cd foam/foam-extend-3.0
source etc/bashrc

( cd wmake/src && make )
cd $WM_THIRD_PARTY_DIR

./AllMake.stage0
./AllMake.stage1
./AllMake.stage2
./AllMake.stage3

cd $WM_PROJECT_DIR
# pick up installed packages
source etc/bashrc

./Allwmake
