#!/usr/bin/python

import pylab
import re
import pickle

forceRegex=r"([0-9.Ee\-+]+)\s+\(+([0-9.Ee\-+]+)\s([0-9.Ee\-+]+)\s([0-9.Ee\-+]+)\)\s\(([0-9.Ee\-+]+)\s([0-9.Ee\-+]+)\s([0-9 .Ee\-+]+)\)+\s\(+([0-9.Ee\-+]+)\s([0-9.Ee\-+]+)\s([0-9.Ee\-+]+)\)\s\(([0-9.Ee\-+]+)\s([0-9.Ee\-+]+)\s([0-9.Ee\-+]+)\)+"
t = []
fpx = []; fpy = []; fpz = []
fvx = []; fvy = []; fvz = []
mpx = []; mpy = []; mpz = []
mvx = []; mvy = []; mvz = []
pipefile=open('forces/0/forces.dat','r')
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

# Calculate total forces
fx = [x + y for x, y in zip(fpx, fvx)]
fy = [x + y for x, y in zip(fpy, fvy)]
fz = [x + y for x, y in zip(fpz, fvz)]

for i in range(len(t)):
    fx[i] = fpx[i] + fvx[i]

for i in range(len(t)):
    fx[i] = fpx[i] + fvx[i]

with open('pressureForces.dat', 'w') as f:
    for f1, f2, f3 in zip(t, fpx, fpy):
        print >> f, f1, f2, f3

with open('viscousForces.dat', 'w') as f:
    for f1, f2, f3 in zip(t, fvx, fvy):
        print >> f, f1, f2, f3

with open('totalForces.dat', 'w') as f:
    for f1, f2, f3 in zip(t, fx, fy):
        print >> f, f1, f2, f3
