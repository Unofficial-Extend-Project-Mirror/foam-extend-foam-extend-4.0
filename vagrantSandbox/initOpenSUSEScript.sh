#! /bin/bash

boxName=$1

echo
echo "Init script for $boxName"
echo

echo
echo "FreeBSD-specific ended. Now doing general stuff"
echo

/vagrant/initGeneralScript.sh

echo
echo "Ended"
