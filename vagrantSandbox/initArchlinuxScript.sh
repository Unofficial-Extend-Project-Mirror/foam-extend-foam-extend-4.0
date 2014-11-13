#! /bin/bash

boxName=$1

echo
echo "Init script for $boxName"
echo

echo
echo "Full update of the system"
echo

pacman -Syu

yaourt --noconfirm -S rpm-org

neededPackages=(gcc-fortran ccache mercurial bison flex git )
bonusPackages=(emacs tcsh)
thirdpartyPackages=(openmpi cmake hwloc)

for p in ${neededPackages[@]}; do
    pacman --noconfirm -S $p
done

for p in ${bonusPackages[@]}; do
    pacman --noconfirm -S $p
done

for p in ${thirdpartyPackages[@]}; do
    pacman --noconfirm -S $p
done

echo
echo "Archlinux-specific ended. Now doing general stuff"
echo

/vagrant/initGeneralScript.sh

echo
echo "Ended"
