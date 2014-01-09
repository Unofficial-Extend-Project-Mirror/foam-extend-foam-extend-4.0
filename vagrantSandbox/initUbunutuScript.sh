#! /bin/bash

boxName=$1

echo
echo "Init script for $boxName"
echo

if [ "$boxName" == "lucid" ]
then
    echo
    echo "Additional Python-Repository"
    echo

    # needed for add-appt-repository
    apt-get -y install python-software-properties

    add-apt-repository ppa:mercurial-ppa/releases
fi

apt-get update -y

echo
echo "Installing additional packages"
echo

apt-get -y install mercurial
apt-get -y install bison
apt-get -y install flex
apt-get -y install g++
apt-get -y install make
#apt-get -y install python-dev
apt-get -y install ccache
apt-get -y install cmake

# test scripts with different shells
apt-get -y install csh
apt-get -y install tcsh
apt-get -y install zsh

# to make the ThirdParty-Stuff work
apt-get -y install rpm

# this is needed for the packaging stuff
echo
echo "Setting for postfix"
echo

# Make sure that default-mta installs
debconf-set-selections <<< "postfix postfix/mailname string vagrant.test.machine.com"
debconf-set-selections <<< "postfix postfix/myhostname string vagrant.test.machine.com"
debconf-set-selections <<< "postfix postfix/main_mailer_type string 'Internet Site'"
debconf-set-selections <<< "postfix postfix/destinations string localhost"

# this workaround doesn't work for lucid
export DEBIAN_FRONTEND=noninteractive

echo
echo "Tools for packaging"
echo

# Needed for packaging
apt-get -y install default-mta
apt-get -y install dpkg-dev
apt-get -y install debhelper devscripts cdbs

# Not needed. Just to keep Bernhard happy
apt-get -y install emacs

echo
echo "Ubuntu-specific ended. Now doing general stuff"
echo

/vagrant/initGeneralScript.sh

echo
echo "Ended"
