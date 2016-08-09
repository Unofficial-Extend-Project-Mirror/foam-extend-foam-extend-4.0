#!/usr/bin/python

import sys
if len(sys.argv) != 2:
        print 'script requires name of log file'
        sys.exit()

logfilename = sys.argv[1]
print 'Reading file', logfilename

import re
UpRegex=r"([A-Z,a-z]*):*.*Solving for Up, Initial residual = \(([0-9.Ee\-+]*)\s([0-9.Ee\-+]*)\s([0-9.Ee\-+]*)\s([0-9.Ee\-+]*)\), Final residual = \(([0-9.Ee\-+]*)\s([0-9.Ee\-+]*)\s([0-9.Ee\-+]*)\s([0-9.Ee\-+]*)\), No Iterations ([0-9]*)"
kRegex=r"([A-Z,a-z]*):*.*Solving for k, Initial residual = ([0-9.Ee\-+]*), Final residual = ([0-9.Ee\-+]*), No Iterations ([0-9]*)"
omegaRegex=r"([A-Z,a-z]*):*.*Solving for omega, Initial residual = ([0-9.Ee\-+]*), Final residual = ([0-9.Ee\-+]*), No Iterations ([0-9]*)"
epsilonRegex=r"([A-Z,a-z]*):*.*Solving for epsilon, Initial residual = ([0-9.Ee\-+]*), Final residual = ([0-9.Ee\-+]*), No Iterations ([0-9]*)"

tUp = []
Ux = []
Uy = []
Uz = []
p = []
iUp = 0

tk = []
k = []
ik = 0

tomega = []
omega = []
iomega = 0

tepsilon = []
epsilon = []
iepsilon = 0

#HJ take name of log file as script argument
pipefile=open(logfilename,'r')
lines = pipefile.readlines()

for line in lines:
        matchUp=re.search(UpRegex,line)
        if matchUp:
                iUp = iUp + 1
                tUp.append(iUp)
                Ux.append(float(matchUp.group(2)))
                Uy.append(float(matchUp.group(3)))
                Uz.append(float(matchUp.group(4)))
                p.append(float(matchUp.group(5)))
        matchk=re.search(kRegex,line)
        if matchk:
                ik = ik + 1
                tk.append(ik)
                k.append(float(matchk.group(2)))
        matchomega=re.search(omegaRegex,line)
        if matchomega:
                iomega = iomega + 1
                tomega.append(iomega)
                omega.append(float(matchomega.group(2)))
        matchepsilon=re.search(epsilonRegex,line)
        if matchepsilon:
                iepsilon = iepsilon + 1
                tepsilon.append(iepsilon)
                epsilon.append(float(matchepsilon.group(2)))

# write clean data file
outfile=open('residual.dat','w')

if iomega > 0:
        for data in zip(tUp,Ux,Uy,Uz,p,k,omega):
                outfile.write(' '.join([str(d) for d in data])+'\n')
elif iepsilon > 0:
        for data in zip(tUp,Ux,Uy,Uz,p,k,epsilon):
                outfile.write(' '.join([str(d) for d in data])+'\n')
elif iUp > 0:
        for data in zip(tUp,Ux,Uy,Uz,p):
                outfile.write(' '.join([str(d) for d in data])+'\n')

outfile.close()

# prepare plot
import pylab
pylab.xlabel('iteration')
pylab.ylabel('residual')
pylab.grid(True)

# plot in log scale
if iUp > 0:
        pylab.semilogy(tUp,Ux,'-',label="Ux")
        pylab.semilogy(tUp,Uy,'-',label="Uy")
        pylab.semilogy(tUp,Uz,'-',label="Uz")
        pylab.semilogy(tUp,p,'-',label="p")

if ik > 0:
        pylab.semilogy(tk,k,'-',label="k")

if iomega > 0:
        pylab.semilogy(tomega,omega,'-',label="omega")

if iepsilon > 0:
        pylab.semilogy(tepsilon,epsilon,'-',label="epsilon")

pylab.legend()
pylab.show()
