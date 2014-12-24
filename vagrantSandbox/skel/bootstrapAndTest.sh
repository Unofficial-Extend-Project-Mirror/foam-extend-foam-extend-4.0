#! /usr/bin/env bash

/home/vagrant/bootstrapFoam.sh

cd foam/foam-extend-3.1
source etc/bashrc

cd testHarness/foam-extend/3.1/runDir

./Allclean

./Allrun_Experimental
