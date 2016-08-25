#! /bin/bash

boxName=$1

echo
echo "Init script for $boxName"
echo

# Otherwise python/mercurial won't install
zypper -n remove patterns-openSUSE-minimal_base-conflicts

# patterns-openSUSE-devel_python

neededPackages=(gcc-c++ mercurial git flex bison make ccache zlib-devel rpm-build cmake)
bonusPackages=(emacs csh tcsh zsh)

for p in ${neededPackages[@]}; do
    zypper -n install $p
done

for p in ${bonusPackages[@]}; do
    zypper -n install $p
done

echo
echo "OpenSUSE-specific ended. Now doing general stuff"
echo

/vagrant/initGeneralScript.sh

echo
echo "Ended"
