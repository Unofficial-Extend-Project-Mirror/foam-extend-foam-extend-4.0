#!/usr/bin/python

import sys
if len(sys.argv) != 2:
        print 'script requires name of log file'
        sys.exit()

logfilename = sys.argv[1]
print 'Reading file', logfilename

import re
UpRegex=r"([A-Z,a-z]*):*.*Solving for Up, Initial residual = \(([0-9.Ee\-+]*)\s([0-9.Ee\-+]*)\s([0-9.Ee\-+]*)\s([0-9.Ee\-+]*)\), Final residual = \(([0-9.Ee\-+]*)\s([0-9.Ee\-+]*)\s([0-9.Ee\-+]*)\s([0-9.Ee\-+]*)\), No Iterations ([0-9]*)"
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
        matchkomega=re.search(komegaRegex,line)
        if matchkomega:
                ikomega = ikomega + 1
                tkomega.append(ikomega)
                k.append(float(matchkomega.group(2)))
                omega.append(float(matchkomega.group(3)))

outfile=open('residual.dat','w')

print 'hits = ', ikomega

#HJ need better way of combining lists
if ikomega > 0:
        for index in range(0,ikomega):
                outfile.write(str(tUp[index])+' '+str(Ux[index])+' '+str(Uy[index])+' '+str(Uz[index])+' '+str(p[index])+' '+str(k[index])+' '+str(omega[index])+'\n')
elif iUp > 0:
        for index in range(0,iUp):
                outfile.write(str(tUp[index])+' '+str(Ux[index])+' '+str(Uy[index])+' '+str(Uz[index])+' '+str(p[index])+'\n')

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

if ikomega > 0:
        pylab.semilogy(tkomega,k,'-',label="k")
        pylab.semilogy(tkomega,omega,'-',label="omega")

pylab.legend()
pylab.show()
