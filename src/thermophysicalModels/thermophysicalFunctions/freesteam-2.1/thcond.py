# -*- coding: utf-8 -*-
# freesteam - IAPWS-IF97 steam tables library
# Copyright (C) 2004-2009  John Pye
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation; either version 2 of
# the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public
# License along with this program; if not, write to the Free
# Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
# Boston, MA 02110-1301, USA.

from math import exp
from freesteam import *

def tc_ptrho(p,T, rho):
	"""
	Revised release on the IAPS Formulation 1985 for the Thermal Conductivity of ordinary water
	IAPWS September 1998
	Page 8

	Converted from Python from the XSteam OpenOffice version.
	"""

	print "rho = %f, T = %f" % (rho, T)

	# ver2.6 Start corrected bug
	if T < 273.15:
		raise RuntimeWarning("T under range")
	elif T < 500 + 273.15:
		if p > 100:
			raise RuntimeWarning("p over range for T<500°C")
	elif T <= 650 + 273.15:
		if p > 70:
			raise RuntimeWarning("p over range for T<650°C")
	elif T <= 800 + 273.15:
		if p > 40:
			raise RuntimeWarning("p over range for T<800°C")
	# ver2.6 End corrected bug

	T = T / 647.26
	rho = rho / 317.7
	lam0 = T ** 0.5 * (0.0102811 + 0.0299621 * T + 0.0156146 * T ** 2 - 0.00422464 * T ** 3)
	lam1 = -0.39707 + 0.400302 * rho + 1.06 * exp(-0.171587 * (rho + 2.39219) ** 2)
	dT = abs(T - 1) + 0.00308976
	Q = 2 + 0.0822994 / dT ** (3. / 5)

	if T >= 1:
		s = 1. / dT
	else:
		s = 10.0932 / dT ** (3. / 5)

	lam2 = (0.0701309 / T ** 10 + 0.011852) * rho ** (9. / 5) * exp(0.642857 * (1 - rho ** (14. / 5))) + 0.00169937 * s * rho ** Q * exp((Q / (1. + Q)) * (1. - rho ** (1. + Q))) - 1.02 * exp(-4.11717 * T ** (3. / 2) - 6.17937 / rho ** 5)

	print "lam0 = %f, lam1 = %f, lam2 = %f" % (lam0, lam1, lam2)

	return lam0 + lam1 + lam2

if __name__=='__main__':
	p = 5
	T = 300
	rho = steam_pT(p * 1e6,T + 273.15).rho
	print "p = %f MPa" % p
	print "T = %f °C" % T
	print "rho = %f kg/m³" % rho
	
	k = tc_ptrho(p, T + 273.15, rho)
	print "k =",k

