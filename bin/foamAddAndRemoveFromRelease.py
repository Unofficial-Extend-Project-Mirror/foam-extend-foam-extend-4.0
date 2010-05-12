#! /usr/bin/python

# debugmode=True
debugmode=False

from os import listdir,path,system
from popen2 import popen4
import sys

def svnCommand(cmd):
    if debugmode:
        print "SVN:",cmd
    else:
        system("svn "+cmd)

def rmEmpty(d):
    if not path.isdir(d):
        return False
    else:
        isEmpty=True
        for f in listdir(d):
            if f==".svn":
                isEmpty=False
            elif not rmEmpty(path.join(d,f)):
                isEmpty=False
        if isEmpty:
            print "Removing ",d,"because it is empty"
            if not debugmode:
                system("rmdir "+d)
        return isEmpty
    
start=sys.argv[1]

rmEmpty(start)

rein,raus=popen4("svn status "+start)
lines=rein.readlines()
rein.close()
raus.close()

modified=0
added=0
removed=0
conflicting=0
replaced=0

for l in lines:
    status=l[0]
    pstatus=l[1]
    name=l[7:-1]
    if status=="?":
        print "Adding",name
        svnCommand("add "+name)
    elif status=="!":
        print "Removing",name
        svnCommand("delete "+name)
    elif status=="M":
        modified+=1
    elif status=="A":
        added+=1
    elif status=="D":
        removed+=1
    elif status=="C":
        conflicting+=1
    elif status=="R":
        replaced+=1
    elif status=="~":
        print "Problem with",name

print
print "Modified files:",modified
print "Added files:",added
print "Removed files:",removed
print "Conflicting files:",conflicting
print "Replaced files:",replaced
print

def checkEmptyDirs(current):
    nrOfContents=0

    for f in listdir(current):
        if f==".svn":
            continue

        pfad=path.join(current,f)

        if path.isdir(pfad):
            if checkEmptyDirs(pfad):
                nrOfContents+=1
        else:
            nrOfContents+=1

    if nrOfContents==0:
        print "Removing",current
        svnCommand("remove "+current)
        return False
    else:
        return True

checkEmptyDirs(start)

