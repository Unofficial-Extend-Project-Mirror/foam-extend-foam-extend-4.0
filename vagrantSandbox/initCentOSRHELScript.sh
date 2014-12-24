#! /bin/bash

boxName=$1

echo
echo "Init script for $boxName"
echo

echo "Install the EPEL-repository for additional software"
rpm -Uvh http://download.fedoraproject.org/pub/epel/6/i386/epel-release-6-8.noarch.rpm

# some of these packages are already installed. But lets be sure

neededPackages=(gcc-c++ gcc-gfortran mercurial git flex bison make ccache rpm-build wget zlib-devel binutils-devel libXt-devel)
bonusPackages=(emacs csh tcsh zsh)

for p in ${neededPackages[@]}; do
    yum install -y $p
done

for p in ${bonusPackages[@]}; do
    yum install -y $p
done

echo "Update mercurial to a more recent version"
rpm -Uhv http://pkgs.repoforge.org/mercurial/mercurial-2.2.2-1.el6.rfx.x86_64.rpm

echo
echo "RHEL/CentOS-specific ended. Now doing general stuff"
echo

/vagrant/initGeneralScript.sh

echo
echo "Ended"
