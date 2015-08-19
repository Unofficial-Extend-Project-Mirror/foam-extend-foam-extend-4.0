#! /usr/bin/env bash

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
        chown -R vagrant:vagrant $target
    fi
done

OFDIR=/home/vagrant/foam/

# make sure that a symbolic link is not erased
if [ ! -e $OFDIR ]; then
    echo "Making directory $OFDIR"
    mkdir -vp $OFDIR
fi

chown -R vagrant:vagrant $OFDIR

# for distros that don't have group vagrant
chown -R vagrant $OFDIR

OFClone=$OFDIR/foam-extend-3.2
OFReference=$OFClone-parent

OFParent=/FOAM-sources

if [ ! -e $OFClone ]
then
    echo
    echo "Cloning the OF-sources"
    echo
    if [ -e "$OFParent/.git" ]
    then
        echo
        echo "Parent is git"
        echo "Cloning. This may take some time"
        echo

        # su -c not correctly working on FreeBSD
        su - vagrant -c "git clone $OFParent $OFClone"

        echo
        echo "Git cloned: TODO: set same branch as parent"
        echo
    elif [ -e "$OFParent/.hg" ]
    then
        echo
        echo "Parent is mercurial. Hello Bernhard"
        echo
#        branchName=`hg branch -R $OFParent`
        idName=`hg id -i -R $OFParent | sed -e "s/\+//"`
        # sed removes + in case of a 'tainted' parent

        echo "Parent is on id $idName"
        echo "Cloning. This may take some time"
        su - vagrant -c "hg clone -u $idName $OFParent $OFClone"
        echo
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
