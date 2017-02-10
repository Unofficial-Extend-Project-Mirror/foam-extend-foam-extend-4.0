#!/usr/bin/python

forcesfilename = 'forces/0/forces.dat'

import sys
if len(sys.argv) != 1:
        print 'script assumes forces file ', forcesfilename
        sys.exit()

import re
forceRegex=r"([0-9.Ee\-+]+)\s+\(+([0-9.Ee\-+]+)\s([0-9.Ee\-+]+)\s([0-9.Ee\-+]+)\)\s\(([0-9.Ee\-+]+)\s([0-9.Ee\-+]+)\s([0-9 .Ee\-+]+)\)+\s\(+([0-9.Ee\-+]+)\s([0-9.Ee\-+]+)\s([0-9.Ee\-+]+)\)\s\(([0-9.Ee\-+]+)\s([0-9.Ee\-+]+)\s([0-9.Ee\-+]+)\)+"
t = []
fpx = []; fpy = []; fpz = []
fvx = []; fvy = []; fvz = []
mpx = []; mpy = []; mpz = []
mvx = []; mvy = []; mvz = []

pipefile=open(forcesfilename,'r')
lines = pipefile.readlines()

for line in lines:
        match=re.search(forceRegex,line)
        if match:
                t.append(float(match.group(1)))
                fpx.append(float(match.group(2)))
                fpy.append(float(match.group(3)))
                fpz.append(float(match.group(4)))
                fvx.append(float(match.group(5)))
                fvy.append(float(match.group(6)))
                fvz.append(float(match.group(7)))
                mpx.append(float(match.group(8)))
                mpy.append(float(match.group(9)))
                mpz.append(float(match.group(10)))
                mvx.append(float(match.group(11)))
                mvy.append(float(match.group(12)))
                mvz.append(float(match.group(13)))


# combine pressure and viscous forces and moments
fx = fpx
fy = fpy
fz = fpz

mx = mpx
my = mpy
mz = mpz

for i in range(1,len(t)):
        fx[i] += fvx[i]
        fy[i] += fvy[i]
        fz[i] += fvz[i]
        mx[i] += mvx[i]
        my[i] += mvy[i]
        mz[i] += mvz[i]

# write clean data file
outForces=open('forces.dat','w')

for data in zip(t,fx,fy,fz):
        outForces.write(' '.join([str(d) for d in data])+'\n')

outForces.close()

outMoments=open('moments.dat','w')

for data in zip(t,mx,my,mz):
        outMoments.write(' '.join([str(d) for d in data])+'\n')

outMoments.close()

# plot forces
import pylab
pylab.xlabel('iteration')
pylab.ylabel('force (N)')
pylab.grid(True)
#
pylab.plot(t,fx,'-',label="fx")
pylab.plot(t,fy,'-',label="fy")
pylab.plot(t,fz,'-',label="fz")

#pylab.plot(t,mx,'-',label="mx")
#pylab.plot(t,my,'-',label="my")
#pylab.plot(t,mz,'-',label="mz")

pylab.legend()
pylab.show()
