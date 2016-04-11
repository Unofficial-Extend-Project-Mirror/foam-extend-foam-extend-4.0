#!/usr/bin/python

import sys
if len(sys.argv) != 2:
        print 'script requires name of log file'
        sys.exit()

logfilename = sys.argv[1]
print 'Reading file', logfilename

import re
UpRegex=r"([A-Z,a-z]*):*.*Solving for Up, Initial residual = \(([0-9.Ee\-+]*)\s([0-9.Ee\-+]*)\s([0-9.Ee\-+]*)\s([0-9.Ee\-+]*)\), Final residual = \(([0-9.Ee\-+]*)\s([0-9.Ee\-+]*)\s([0-9.Ee\-+]*)\s([0-9.Ee\-+]*)\), No Iterations ([0-9]*)"
kepsilonRegex=r"([A-Z,a-z]*):*.*Solving for kEpsilon, Initial residual = \(([0-9.Ee\-+]*)\s([0-9.Ee\-+]*)\), Final residual = \(([0-9.Ee\-+]*)\s([0-9.Ee\-+]*)\), No Iterations ([0-9]*)"
komegaRegex=r"([A-Z,a-z]*):*.*Solving for kOmega, Initial residual = \(([0-9.Ee\-+]*)\s([0-9.Ee\-+]*)\), Final residual = \(([0-9.Ee\-+]*)\s([0-9.Ee\-+]*)\), No Iterations ([0-9]*)"

tUp = []
Ux = []
Uy = []
Uz = []
p = []
iUp = 0

tkomega = []
k = []
omega = []
ikomega = 0

tkepsilon = []
k = []
epsilon = []
ikepsilon = 0

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
        matchkepsilon=re.search(kepsilonRegex,line)
        if matchkepsilon:
                ikepsilon = ikepsilon + 1
                tkepsilon.append(ikepsilon)
                k.append(float(matchkepsilon.group(2)))
                epsilon.append(float(matchkepsilon.group(3)))
        matchkomega=re.search(komegaRegex,line)
        if matchkomega:
                ikomega = ikomega + 1
                tkomega.append(ikomega)
                k.append(float(matchkomega.group(2)))
                omega.append(float(matchkomega.group(3)))

outfile=open('residual.dat','w')

if ikomega > 0:
        for data in zip(tUp,Ux,Uy,Uz,p,k,omega):
                outfile.write(' '.join([str(d) for d in data])+'\n')
elif ikepsilon > 0:
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

if ikepsilon > 0:
        pylab.semilogy(tkepsilon,k,'-',label="k")
        pylab.semilogy(tkepsilon,epsilon,'-',label="epsilon")

if ikomega > 0:
        pylab.semilogy(tkomega,k,'-',label="k")
        pylab.semilogy(tkomega,omega,'-',label="omega")

pylab.legend()
pylab.show()
