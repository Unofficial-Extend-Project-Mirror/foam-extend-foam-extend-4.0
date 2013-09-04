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

    apt-get update -y
fi

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
echo "Copying stuff from skeleton"
echo
for f in $(ls -A /vagrant/skel)
do
    target="/home/vagrant/$f"
    if [ -e $target ]
    then
	echo "$target already there"
    else
	echo "Copying $target from skeleton"
	cp -r "/vagrant/skel/$f" $target
    fi
done

OFDIR=/home/vagrant/OpenFOAM/

mkdir -vp $OFDIR

OFClone=$OFDIR/OpenFOAM-1.6-ext
OFReference=$OFDIR/OpenFOAM-1.6-ext-parent

OFParent=/OpenFOAM-sources

if [ ! -e $OFClone ]
then
    echo
    echo "Cloning the OF-sources"
    echo
    if [ -e "$OFParent/.git" ]
    then
	echo
	echo "Parent is git"
	echo
	su vagrant - -c "git clone $OFParent $OFClone"
	echo
	echo "Git cloned: TODO: set same branch as parent"
	echo
    elif [ -e "$OFParent/.hg" ]
    then
	echo
	echo "Parent is mercurial. Hello Bernhard"
	echo
	su vagrant - -c "hg clone $OFParent $OFClone"
	branchName=`hg branch -R $OFParent`
	echo
	echo "Parent is on branch $branchName"
	su vagrant - -c "hg update -R $OFClone $branchName"
    else
	echo
	echo "Problem. Parent $OFParent is neither git nor mercurial"
	echo
    fi
else
    echo "Repository $OFClone already there. No cloning"
fi

if [ ! -e $OFReference ]
then
    echo
    echo "Linking $OFReference to $OFParent"
    echo
    ln -s $OFParent $OFReference
else
    echo
    echo "Link $OFReference already there"
    echo
fi

chown -R vagrant:vagrant $OFDIR

echo
echo "Current ccache:"
export CCACHE_DIR=/vagrant/ccache4vm; ccache --show-stats

echo
echo "Ended"
