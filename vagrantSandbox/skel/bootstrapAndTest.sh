#! /usr/bin/env bash

/home/vagrant/bootstrapFoam.sh

cd foam/foam-extend-3.2
source etc/bashrc

cd testHarness/foam-extend/3.2/runDir

./Allclean

./Allrun_Experimental
