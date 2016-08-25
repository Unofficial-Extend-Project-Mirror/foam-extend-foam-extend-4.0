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

neededPackages=(g++ bison flex mercurial git make ccache cmake rpm zlib1g-dev libiberty-dev)
bonusPackages=(emacs csh tcsh zsh)

for p in ${neededPackages[@]}; do
    apt-get -y install $p
done

for p in ${bonusPackages[@]}; do
    apt-get -y install $p
done

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

packagingPackages=(default-mta dpkg-dev debhelper devscripts cdbs)
for p in ${packagingPackages[@]}; do
    apt-get -y install $p
done

echo
echo "Ubuntu-specific ended. Now doing general stuff"
echo

/vagrant/initGeneralScript.sh

echo
echo "Ended"
