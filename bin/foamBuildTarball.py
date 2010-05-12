#! /usr/bin/python

import sys,re
from os import path
from subprocess import Popen,PIPE,call
import tarfile 

if len(sys.argv)!=2:
    print "Error: SVN-Url is needed"
    sys.exit(-1)

url=sys.argv[1]

name=path.basename(url[:-1])

p=Popen(["svn","info",url],stdin=PIPE, stdout=PIPE, close_fds=True)

(child_stdout, child_stdin) = (p.stdout, p.stdin)

revision=-1

for l in child_stdout.readlines():
    m=re.compile("Last Changed Rev: (.+)").match(l)
    if m!=None:
        revision=int(m.group(1))

if revision<0:
    print "Invalid URL or stuff"
    sys.exit(-1)
    
fullname="%s.r%d" % (name,revision)
l
print "Generating",fullname

retcode=call(["svn","export",url,fullname])
if retcode!=0:
    print "Problem. Returncode",retcode
    sys.exit(-1)

print "Tarring ...."
tar=tarfile.open(fullname+".tgz","w:gz")
tar.add(fullname,arcname=name)
tar.close()
print "Removing directory"
retcode=call(["rm","-rf",fullname])
print "Finished"
