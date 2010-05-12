#!/usr/bin/python

# this script adds a set of SVN-properties to files and directories under
# a directory that is specified on the command line

from popen2 import popen2
import sys
import string
import glob
from os import path,listdir

svnCommand="svn"
isSVK=False

def runSvn(cmd):
    raus,rein=popen2(svnCommand+" "+cmd)
    result=raus.readlines()
    rein.close()
    raus.close()
    return result

def getProperty(fName,property):
    raw=runSvn("propget %s %s" % (property,fName))
    return string.join(raw)

def setProperty(fName,property,value):
    runSvn("propset %s \"%s\" %s" % (property,value,fName))

def addToListProperty(fName,property,value):
    tmp=getProperty(fName,property)
    lst=map(string.strip,string.split(tmp))
    if not value in lst:
        lst.append(value)
    else:
        return False
    val=string.join(lst,"\n")
    setProperty(fName,property,val)
    return True

def addKeyword(fName,keyword):
    return addToListProperty(fName,"svn:keywords",keyword)
    
def addIgnore(fName,keyword):
    return addToListProperty(fName,"svn:ignore",keyword)
    
def recursivlyDoToFiles(directory,fileFilter,function,isDir=False,testSvn=True):
    if testSvn and not isSVK:
        if not path.exists(path.join(directory,".svn")):
            return
        
    for f in glob.glob(path.join(directory,fileFilter)):
        if not path.isfile(f) and not path.isdir(f):
            continue
        
        if (isDir and path.isfile(f)) or (not isDir and path.isdir(f)):
            continue
        
        if isDir and testSvn and not isSVK:
            if not path.exists(path.join(f,".svn")):
                continue
                
        if function(f):
            print "....",f
                
    for f in listdir(directory):
        if f not in [".svn","lnInclude"]:
            tmp=path.join(directory,f)
            if path.isdir(tmp):
                recursivlyDoToFiles(tmp,fileFilter,function,isDir=isDir,testSvn=testSvn)

if not path.exists(path.join(sys.argv[1],".svn")):
    svnCommand="svk"
    isSVK=True
    
print "\nAdding Id-keyword to Python-files"
recursivlyDoToFiles(sys.argv[1],"*.py",lambda x:addKeyword(x,"Id"))

print "\nAdding Id-keyword to C++-files"
recursivlyDoToFiles(sys.argv[1],"*.C",lambda x:addKeyword(x,"Id"))

print "\nAdding Id-keyword to C++-headers"
recursivlyDoToFiles(sys.argv[1],"*.H",lambda x:addKeyword(x,"Id"))

print "\nAdding *Opt to ignore-list for Make-directories"
recursivlyDoToFiles(sys.argv[1],"Make",lambda x:addIgnore(x,"*Opt"),isDir=True)

print "\nAdding *Debug to ignore-list for Make-directories"
recursivlyDoToFiles(sys.argv[1],"Make",lambda x:addIgnore(x,"*Debug"),isDir=True)

print "\nAdding lnInclude to ignore-list for all directories"
recursivlyDoToFiles(sys.argv[1],"*",lambda x:addIgnore(x,"lnInclude"),isDir=True)

print "\nAdding *.dep to ignore-list for all directories"
recursivlyDoToFiles(sys.argv[1],"*",lambda x:addIgnore(x,"*.dep"),isDir=True)
