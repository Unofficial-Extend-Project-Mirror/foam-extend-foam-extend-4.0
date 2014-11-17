#! /usr/bin/env bash

boxName=$1

echo
echo "Init script for $boxName"
echo

pkg install -fy virtualbox-ose-additions

neededPackages=(mercurial git flex bison ccache rpm4 wget)
bonusPackages=(emacs24  zsh)

for p in ${neededPackages[@]}; do
    pkg install -y $p
done

for p in ${bonusPackages[@]}; do
    pkg install -y $p
done

echo
echo "Upgrading all to get working packages"
echo
pkg upgrade -y

echo
echo "FreeBSD-specific ended. Now doing general stuff"
echo

/vagrant/initGeneralScript.sh

echo
echo "Ended"
